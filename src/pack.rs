//! Lipid coordinate relaxation: a centroidal-Voronoi (Lloyd) pass that drives the
//! lateral density toward uniform, followed by a repulsive de-clash so the start is
//! overlap-free.
//!
//! The build seeds each leaflet/shell with a **random** placement; the routines
//! here turn that into a clean, uniform, clash-free starting structure. The Voronoi
//! relaxation lives in [`crate::voronoi`] (per-site clipping); this module supplies
//! the soft-sphere de-clash and orchestrates the two:
//!
//! - [`equalize_plane`] for a flat bilayer: a per-leaflet planar Lloyd pass then a
//!   de-clash, in a periodic box.
//! - [`equalize_on_sphere`] for a vesicle: a per-shell spherical Lloyd pass then a
//!   de-clash.
//! - [`equalize_on_cylinder`] for a tubule: a per-shell cylindrical Lloyd pass
//!   (isometrically unrolled to a periodic plane) then a cylindrical de-clash.
//!
//! Beads are treated as soft spheres: an overlapping pair pushes apart by a fraction
//! of the overlap each step (a damped Jacobi relaxation), capped for stability, over
//! a flat CSR cell list parallelised per lipid. A whole lipid moves rigidly
//! (laterally in the plane, with z fixed so it stays in its leaflet, tangentially
//! over the sphere at fixed radius, or tangentially + axially on the cylinder at
//! fixed radius), so the relaxation never distorts a lipid or moves it out of its
//! layer.

use std::f64::consts::TAU;
use std::ops::Range;

use rayon::prelude::*;

use crate::voronoi::{
    PlaneCvtWorkspace, SphereCvtWorkspace, extract_plane_sites, extract_sphere_dirs,
    rotate_lipid_from_to, rotate_with_sin_cos, wrap,
};

/// Fraction of each pair's overlap removed per step (under-relaxed for stability).
const STEP_FRACTION: f64 = 0.5;
/// Per-step lipid displacement cap, as a fraction of the largest bead radius.
const MAX_STEP_FRACTION: f64 = 0.3;
/// Stop a de-clash once the largest overlap stops shrinking by at least this much
/// (nm) for [`STALL_LIMIT`] straight steps: a stuck contact the slide cannot remove.
const MIN_IMPROVEMENT: f64 = 1e-4;
const STALL_LIMIT: usize = 25;

/// Lloyd (CVT) round cap before the final de-clash. A few dozen rounds give a
/// uniform-enough start; the equilibration MD refines it further. Convergence
/// usually stops the loop before the cap.
const CVT_ROUNDS: usize = 60;
/// CVT convergence: stop once the largest lipid move in a round falls below this
/// fraction of the mean lipid spacing.
const CVT_CONVERGED_FRACTION: f64 = 1e-2;

/// Minimum chunk size for Rayon parallel lipid sweeps.
const PAR_MIN_LEN: usize = 64;

/// Packed bead data stored contiguously in cell order for cache-friendly pair scans.
#[derive(Clone, Copy, Default)]
struct PackBead {
    pos: [f64; 3],
    r: f64,
    lipid: u32,
}

/// Shared convergence / stall detector for both planar and spherical de-clash loops.
struct StallTracker {
    best: f64,
    stalled: usize,
}

impl StallTracker {
    fn new() -> Self {
        Self {
            best: f64::MAX,
            stalled: 0,
        }
    }

    fn should_stop(&mut self, max_overlap: f64, tolerance: f64) -> bool {
        if max_overlap < tolerance {
            return true;
        }
        if self.best - max_overlap < MIN_IMPROVEMENT {
            self.stalled += 1;
            if self.stalled >= STALL_LIMIT {
                return true;
            }
        } else {
            self.stalled = 0;
        }
        self.best = self.best.min(max_overlap);
        false
    }
}

/// Build the per-bead owning lipid index (`u32::MAX` for unowned beads).
fn bead_owners(n_beads: usize, lipids: &[Range<usize>]) -> Vec<u32> {
    let mut owner = vec![u32::MAX; n_beads];
    for (lipid, range) in lipids.iter().enumerate() {
        owner[range.clone()].fill(lipid as u32);
    }
    owner
}

/// Prefix sum `counts` into `offsets`, resetting `counts` to serve as CSR write cursors.
fn prefix_sum(counts: &mut [u32], offsets: &mut [u32]) {
    let mut acc = 0u32;
    for (c, o) in counts.iter_mut().zip(offsets.iter_mut()) {
        *o = acc;
        let len = *c;
        *c = acc;
        acc += len;
    }
    offsets[counts.len()] = acc;
}

/// Wrap cell index `g + d` (for `d ∈ {-1, 0, 1}`) modulo `n` without division.
#[inline]
fn wrap_cell(g: u32, d: isize, n: usize) -> usize {
    match d {
        -1 => {
            if g == 0 {
                n - 1
            } else {
                g as usize - 1
            }
        }
        0 => g as usize,
        1 => {
            let next = g as usize + 1;
            if next == n { 0 } else { next }
        }
        _ => (g as isize + d).rem_euclid(n as isize) as usize,
    }
}

// ============================ Planar ============================

/// A fixed spherical obstacle (e.g. a coarse-grained protein bead plus
/// clearance margin) that lipids must flow around during CVT and de-clash.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Obstacle {
    /// Cartesian position (nm) in the same coordinate frame as `positions`.
    pub pos: [f64; 3],
    /// Effective excluded radius (nm).
    pub radius: f64,
}

/// Minimum spacing (nm) between projected pinned CVT obstacle sites so dense
/// protein cores form a clean boundary without redundant duplicates.
const PINNED_SITE_SPACING: f64 = 0.35;

/// Deduplicate 2D periodic pinned sites closer than `min_dist` (nm).
///
/// Uses an O(N) periodic cell grid with an intrusive linked list for large
/// inputs, preserving exact greedy insertion order.
fn dedup_plane_sites(candidates: Vec<[f64; 2]>, box_xy: [f64; 2], min_dist: f64) -> Vec<[f64; 2]> {
    let min_d2 = min_dist * min_dist;
    let half_x = 0.5 * box_xy[0];
    let half_y = 0.5 * box_xy[1];
    let periodic_close = |p: [f64; 2], q: [f64; 2]| -> bool {
        let mut dx = p[0] - q[0];
        if dx > half_x {
            dx -= box_xy[0];
        } else if dx < -half_x {
            dx += box_xy[0];
        }
        let mut dy = p[1] - q[1];
        if dy > half_y {
            dy -= box_xy[1];
        } else if dy < -half_y {
            dy += box_xy[1];
        }
        dx * dx + dy * dy < min_d2
    };

    let mut kept: Vec<[f64; 2]> = Vec::new();
    let nx = ((box_xy[0] / min_dist).floor() as usize).clamp(1, 512);
    let ny = ((box_xy[1] / min_dist).floor() as usize).clamp(1, 512);
    if candidates.len() <= 64 || nx.saturating_mul(ny) > 262_144 {
        for p in candidates {
            if !kept.iter().any(|&q| periodic_close(p, q)) {
                kept.push(p);
            }
        }
        return kept;
    }

    let inv_wx = nx as f64 / box_xy[0];
    let inv_wy = ny as f64 / box_xy[1];
    let dxs: &[isize] = match nx {
        1 => &[0],
        2 => &[0, 1],
        _ => &[-1, 0, 1],
    };
    let dys: &[isize] = match ny {
        1 => &[0],
        2 => &[0, 1],
        _ => &[-1, 0, 1],
    };
    let mut head = vec![u32::MAX; nx * ny];
    let mut next: Vec<u32> = Vec::with_capacity(candidates.len());

    for p in candidates {
        let cx = ((p[0] * inv_wx) as usize).min(nx - 1);
        let cy = ((p[1] * inv_wy) as usize).min(ny - 1);
        let mut too_close = false;
        'search: for &dy in dys {
            let ny_i = (cy as isize + dy).rem_euclid(ny as isize) as usize;
            let row = ny_i * nx;
            for &dx in dxs {
                let nx_i = (cx as isize + dx).rem_euclid(nx as isize) as usize;
                let mut cur = head[row + nx_i];
                while cur != u32::MAX {
                    if periodic_close(p, kept[cur as usize]) {
                        too_close = true;
                        break 'search;
                    }
                    cur = next[cur as usize];
                }
            }
        }
        if !too_close {
            let idx = kept.len() as u32;
            kept.push(p);
            let cell = cy * nx + cx;
            next.push(head[cell]);
            head[cell] = idx;
        }
    }
    kept
}

/// Deduplicate unit sphere directions closer than `min_chord` in 3D chord distance.
///
/// Uses an O(N) 3D cell grid over `[-1, 1]³` for large inputs.
fn dedup_sphere_dirs(candidates: Vec<[f64; 3]>, min_chord: f64) -> Vec<[f64; 3]> {
    let min_d2 = min_chord * min_chord;
    let chord_close = |d: [f64; 3], q: [f64; 3]| -> bool {
        let dx = d[0] - q[0];
        let dy = d[1] - q[1];
        let dz = d[2] - q[2];
        dx * dx + dy * dy + dz * dz < min_d2
    };

    let mut kept: Vec<[f64; 3]> = Vec::new();
    if candidates.len() <= 64 || min_chord <= 0.0 {
        for d in candidates {
            if !kept.iter().any(|&q| chord_close(d, q)) {
                kept.push(d);
            }
        }
        return kept;
    }

    // Cell width `2 / n >= min_chord`, so the 3×3×3 cell neighborhood is exact.
    let n = ((2.0 / min_chord).floor() as usize).clamp(1, 32);
    let inv_w = 0.5 * n as f64;
    let mut head = vec![u32::MAX; n * n * n];
    let mut next: Vec<u32> = Vec::with_capacity(candidates.len());

    for d in candidates {
        let cx = (((d[0] + 1.0) * inv_w) as usize).min(n - 1);
        let cy = (((d[1] + 1.0) * inv_w) as usize).min(n - 1);
        let cz = (((d[2] + 1.0) * inv_w) as usize).min(n - 1);
        let mut too_close = false;
        'search: for iz in cz.saturating_sub(1)..=(cz + 1).min(n - 1) {
            let z_base = iz * n;
            for iy in cy.saturating_sub(1)..=(cy + 1).min(n - 1) {
                let row = (z_base + iy) * n;
                for ix in cx.saturating_sub(1)..=(cx + 1).min(n - 1) {
                    let mut cur = head[row + ix];
                    while cur != u32::MAX {
                        if chord_close(d, kept[cur as usize]) {
                            too_close = true;
                            break 'search;
                        }
                        cur = next[cur as usize];
                    }
                }
            }
        }
        if !too_close {
            let idx = kept.len() as u32;
            kept.push(d);
            let cell = (cz * n + cy) * n + cx;
            next.push(head[cell]);
            head[cell] = idx;
        }
    }
    kept
}

/// An observer called with the current bead positions after each relaxation step
/// (each Lloyd round, then each de-clash iteration), so the binary can record
/// frames for the `--gif` animation. `None`, the usual case, adds no work.
pub type StepObserver<'a> = &'a mut dyn FnMut(&[[f64; 3]]);

/// Centroidal-Voronoi relax a flat bilayer to a uniform density, then de-clash.
///
/// `leaflets` lists each leaflet's lipid ranges (relaxed as an independent periodic
/// CVT at its own density); `positions` holds every lipid bead (mutated in place,
/// x/y only); `radii` is the de-clash radius per bead. Returns the largest remaining
/// bead overlap (nm) after the de-clash; near zero is a clash-free start.
pub fn equalize_plane(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    leaflets: &[&[Range<usize>]],
    box_xy: [f64; 2],
    iterations: usize,
    tolerance: f64,
    observer: Option<StepObserver>,
) -> f64 {
    equalize_plane_with_obstacles(
        positions,
        radii,
        leaflets,
        &[],
        box_xy,
        iterations,
        tolerance,
        observer,
    )
}

/// Like [`equalize_plane`], but treats `obstacles` as fixed spheres that clip
/// nearby lipids' Voronoi cells and repel overlapping lipid beads during de-clash.
#[allow(clippy::too_many_arguments)]
pub fn equalize_plane_with_obstacles(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    leaflets: &[&[Range<usize>]],
    obstacles: &[Obstacle],
    box_xy: [f64; 2],
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let all: Vec<Range<usize>> = leaflets.iter().flat_map(|l| l.iter().cloned()).collect();
    let total: usize = all.len();
    if total == 0 || box_xy[0] <= 0.0 || box_xy[1] <= 0.0 {
        return 0.0;
    }
    let converged = CVT_CONVERGED_FRACTION * (box_xy[0] * box_xy[1] / total as f64).sqrt();
    let mut cvt = PlaneCvtWorkspace::new();

    let max_lipid_r = radii.iter().copied().fold(0.0_f64, f64::max);
    // Run CVT directly on per-lipid 2D sites, appending deduplicated pinned
    // obstacle sites that overlap each leaflet's z-slab.
    let init_sites: Vec<Vec<[f64; 2]>> = leaflets
        .iter()
        .map(|l| extract_plane_sites(positions, l, box_xy))
        .collect();
    let mut cur_sites: Vec<(usize, Vec<[f64; 2]>)> = leaflets
        .iter()
        .zip(&init_sites)
        .map(|(l, sites)| {
            let n_mov = sites.len();
            let mut combined = sites.clone();
            if !obstacles.is_empty() && !l.is_empty() {
                let (mut z_lo, mut z_hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for range in *l {
                    for k in range.clone() {
                        z_lo = z_lo.min(positions[k][2]);
                        z_hi = z_hi.max(positions[k][2]);
                    }
                }
                let raw: Vec<[f64; 2]> = obstacles
                    .iter()
                    .filter(|o| {
                        let reach = o.radius + max_lipid_r;
                        o.pos[2] + reach >= z_lo && o.pos[2] - reach <= z_hi
                    })
                    .map(|o| [wrap(o.pos[0], box_xy[0]), wrap(o.pos[1], box_xy[1])])
                    .collect();
                combined.extend(dedup_plane_sites(raw, box_xy, PINNED_SITE_SPACING));
            }
            (n_mov, combined)
        })
        .collect();
    let record_each = observer.is_some();

    for _ in 0..CVT_ROUNDS {
        let mut moved = 0.0_f64;
        for (leaflet, (n_mov, sites)) in leaflets.iter().zip(&mut cur_sites) {
            moved = moved.max(cvt.step_sites_pinned(sites, *n_mov, box_xy));
            if record_each {
                for (range, &[dx, dy]) in leaflet.iter().zip(cvt.last_displacements()) {
                    for k in range.clone() {
                        positions[k][0] = wrap(positions[k][0] + dx, box_xy[0]);
                        positions[k][1] = wrap(positions[k][1] + dy, box_xy[1]);
                    }
                }
            }
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
        if moved < converged {
            break;
        }
    }
    if !record_each {
        for (leaflet, (init, (_, cur))) in leaflets.iter().zip(init_sites.iter().zip(&cur_sites)) {
            for (range, (&s0, &s1)) in leaflet.iter().zip(init.iter().zip(cur)) {
                let dx = s1[0] - s0[0];
                let dy = s1[1] - s0[1];
                for k in range.clone() {
                    positions[k][0] = wrap(positions[k][0] + dx, box_xy[0]);
                    positions[k][1] = wrap(positions[k][1] + dy, box_xy[1]);
                }
            }
        }
    }
    relax(
        positions, radii, &all, obstacles, box_xy, iterations, tolerance, observer,
    )
}

/// De-clash `lipids` laterally so no bead overlaps another lipid's bead or any
/// fixed `obstacles`, in a periodic box `box_xy` (nm). Each lipid translates
/// rigidly in x/y; z is fixed. Stops when the largest overlap drops below
/// `tolerance` (nm) or after `iterations`. Returns the largest remaining overlap.
#[allow(clippy::too_many_arguments)]
fn relax(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    lipids: &[Range<usize>],
    obstacles: &[Obstacle],
    box_xy: [f64; 2],
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let max_radius = radii
        .iter()
        .copied()
        .chain(obstacles.iter().map(|o| o.radius))
        .fold(0.0_f64, f64::max);
    if max_radius <= 0.0 || lipids.is_empty() {
        return 0.0;
    }
    let n_beads = positions.len();
    let n_obs = obstacles.len();
    let n_total = n_beads + n_obs;
    let owner = bead_owners(n_beads, lipids);

    // Cells at least one bead diameter wide, so any overlapping pair is in adjacent
    // cells. The grid covers the periodic x/y box.
    let (bx, by) = (box_xy[0], box_xy[1]);
    let (half_bx, half_by) = (0.5 * bx, 0.5 * by);
    let cell = 2.0 * max_radius;
    let nx = ((bx / cell).floor() as usize).max(1);
    let ny = ((by / cell).floor() as usize).max(1);
    let inv_wx = nx as f64 / bx;
    let inv_wy = ny as f64 / by;
    let ncells = nx * ny;

    // Neighbour-cell offsets per axis without double-counting when nx or ny <= 2.
    let dxs: &[isize] = match nx {
        1 => &[0],
        2 => &[0, 1],
        _ => &[-1, 0, 1],
    };
    let dys: &[isize] = match ny {
        1 => &[0],
        2 => &[0, 1],
        _ => &[-1, 0, 1],
    };
    let max_step = MAX_STEP_FRACTION * max_radius;

    // Pre-wrap fixed obstacles once.
    let obs_wrapped: Vec<([f64; 3], u32, u32, f64)> = obstacles
        .iter()
        .map(|o| {
            let wx = wrap(o.pos[0], bx);
            let wy = wrap(o.pos[1], by);
            let gx = ((wx * inv_wx) as usize).min(nx - 1) as u32;
            let gy = ((wy * inv_wy) as usize).min(ny - 1) as u32;
            ([wx, wy, o.pos[2]], gx, gy, o.radius)
        })
        .collect();

    // Reusable CSR cell-list buffers across iterations (zero hot-loop allocs).
    let mut wrapped = vec![[0.0_f64; 3]; n_beads];
    let mut bead_cell = vec![(0u32, 0u32); n_beads];
    let mut counts = vec![0u32; ncells];
    let mut offsets = vec![0u32; ncells + 1];
    let mut sorted = vec![PackBead::default(); n_total];
    let mut force = vec![[0.0_f64; 2]; lipids.len()];

    let mut max_overlap = 0.0_f64;
    let mut tracker = StallTracker::new();
    for _ in 0..iterations {
        counts.fill(0);
        for (i, &p) in positions.iter().enumerate() {
            let wx = wrap(p[0], bx);
            let wy = wrap(p[1], by);
            let gx = ((wx * inv_wx) as usize).min(nx - 1);
            let gy = ((wy * inv_wy) as usize).min(ny - 1);
            wrapped[i] = [wx, wy, p[2]];
            bead_cell[i] = (gx as u32, gy as u32);
            counts[gy * nx + gx] += 1;
        }
        for &(_, gx, gy, _) in &obs_wrapped {
            counts[gy as usize * nx + gx as usize] += 1;
        }
        prefix_sum(&mut counts, &mut offsets);
        for i in 0..n_beads {
            let (gx, gy) = bead_cell[i];
            let c = gy as usize * nx + gx as usize;
            let slot = counts[c] as usize;
            counts[c] += 1;
            sorted[slot] = PackBead {
                pos: wrapped[i],
                r: radii[i],
                lipid: owner[i],
            };
        }
        for &(pos, gx, gy, r) in &obs_wrapped {
            let c = gy as usize * nx + gx as usize;
            let slot = counts[c] as usize;
            counts[c] += 1;
            sorted[slot] = PackBead {
                pos,
                r,
                lipid: u32::MAX,
            };
        }

        let wrapped_ref = &wrapped;
        let bead_cell_ref = &bead_cell;
        let offsets_ref = &offsets;
        let sorted_ref = &sorted;

        max_overlap = lipids
            .par_iter()
            .enumerate()
            .zip(force.par_iter_mut())
            .with_min_len(PAR_MIN_LEN)
            .map(|((li, range), f_out)| {
                let li_u32 = li as u32;
                let mut fx = 0.0_f64;
                let mut fy = 0.0_f64;
                let mut local_max = 0.0_f64;
                for i in range.clone() {
                    let pi = wrapped_ref[i];
                    let ri = radii[i];
                    let (gx, gy) = bead_cell_ref[i];
                    for &dyc in dys {
                        let row = wrap_cell(gy, dyc, ny) * nx;
                        for &dxc in dxs {
                            let c = row + wrap_cell(gx, dxc, nx);
                            let start = offsets_ref[c] as usize;
                            let end = offsets_ref[c + 1] as usize;
                            for bj in &sorted_ref[start..end] {
                                if bj.lipid == li_u32 {
                                    continue;
                                }
                                let sigma = ri + bj.r;
                                let sigma2 = sigma * sigma;
                                let dz = pi[2] - bj.pos[2];
                                let dz2 = dz * dz;
                                if dz2 >= sigma2 {
                                    continue;
                                }
                                let mut dx = pi[0] - bj.pos[0];
                                if dx > half_bx {
                                    dx -= bx;
                                } else if dx < -half_bx {
                                    dx += bx;
                                }
                                let mut dy = pi[1] - bj.pos[1];
                                if dy > half_by {
                                    dy -= by;
                                } else if dy < -half_by {
                                    dy += by;
                                }
                                let d2 = dx * dx + dy * dy + dz2;
                                if d2 < sigma2 && d2 > 1e-12 {
                                    let d = d2.sqrt();
                                    let overlap = sigma - d;
                                    local_max = local_max.max(overlap);
                                    let scale = overlap / d;
                                    fx += dx * scale;
                                    fy += dy * scale;
                                }
                            }
                        }
                    }
                }
                *f_out = [fx, fy];
                local_max
            })
            .reduce(|| 0.0_f64, f64::max);

        if tracker.should_stop(max_overlap, tolerance) {
            break;
        }

        for (range, &[fx, fy]) in lipids.iter().zip(&force) {
            if fx == 0.0 && fy == 0.0 {
                continue;
            }
            let mut sx = STEP_FRACTION * fx;
            let mut sy = STEP_FRACTION * fy;
            let mag = (sx * sx + sy * sy).sqrt();
            if mag > max_step {
                let scale = max_step / mag;
                sx *= scale;
                sy *= scale;
            }
            for i in range.clone() {
                positions[i][0] += sx;
                positions[i][1] += sy;
            }
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
    }
    max_overlap
}

// ============================ Spherical ============================

/// Centroidal-Voronoi relax a vesicle to a uniform areal density, then de-clash.
///
/// `shells` pairs each shell's lipid ranges with its radius (nm); each is relaxed as
/// an independent spherical CVT. `positions` holds every lipid bead (mutated in
/// place); `radii` is the de-clash radius per bead. Returns the largest remaining
/// bead overlap (nm).
pub fn equalize_on_sphere(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    shells: &[(&[Range<usize>], f64)],
    center: [f64; 3],
    iterations: usize,
    tolerance: f64,
    observer: Option<StepObserver>,
) -> f64 {
    equalize_on_sphere_with_obstacles(
        positions,
        radii,
        shells,
        &[],
        center,
        iterations,
        tolerance,
        observer,
    )
}

/// Like [`equalize_on_sphere`], but treats `obstacles` as fixed spheres that
/// clip nearby spherical Voronoi cells and repel overlapping lipid beads.
#[allow(clippy::too_many_arguments)]
pub fn equalize_on_sphere_with_obstacles(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    shells: &[(&[Range<usize>], f64)],
    obstacles: &[Obstacle],
    center: [f64; 3],
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let all: Vec<Range<usize>> = shells.iter().flat_map(|(s, _)| s.iter().cloned()).collect();
    if all.is_empty() {
        return 0.0;
    }
    let converged = shells
        .iter()
        .filter(|(s, _)| !s.is_empty())
        .map(|(s, r)| {
            CVT_CONVERGED_FRACTION * r * (4.0 * std::f64::consts::PI / s.len() as f64).sqrt()
        })
        .fold(f64::MAX, f64::min);
    let mut cvt = SphereCvtWorkspace::new();

    let max_lipid_r = radii.iter().copied().fold(0.0_f64, f64::max);
    let init_dirs: Vec<Vec<[f64; 3]>> = shells
        .iter()
        .map(|(s, _)| extract_sphere_dirs(positions, s, center))
        .collect();
    let mut cur_dirs: Vec<(usize, Vec<[f64; 3]>)> = shells
        .iter()
        .zip(&init_dirs)
        .map(|((shell, radius), dirs)| {
            let n_mov = dirs.len();
            let mut combined = dirs.clone();
            if !obstacles.is_empty() && !shell.is_empty() && *radius > 0.0 {
                let (mut r_lo, mut r_hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for range in *shell {
                    for k in range.clone() {
                        let d = [
                            positions[k][0] - center[0],
                            positions[k][1] - center[1],
                            positions[k][2] - center[2],
                        ];
                        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        r_lo = r_lo.min(r);
                        r_hi = r_hi.max(r);
                    }
                }
                let raw: Vec<[f64; 3]> = obstacles
                    .iter()
                    .filter_map(|o| {
                        let d = [
                            o.pos[0] - center[0],
                            o.pos[1] - center[1],
                            o.pos[2] - center[2],
                        ];
                        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        let reach = o.radius + max_lipid_r;
                        (r > 1e-6 && r + reach >= r_lo && r - reach <= r_hi).then_some([
                            d[0] / r,
                            d[1] / r,
                            d[2] / r,
                        ])
                    })
                    .collect();
                let min_chord = (PINNED_SITE_SPACING / *radius).min(0.5);
                combined.extend(dedup_sphere_dirs(raw, min_chord));
            }
            (n_mov, combined)
        })
        .collect();
    let record_each = observer.is_some();

    for _ in 0..CVT_ROUNDS {
        let mut moved = 0.0_f64;
        for ((shell, radius), (n_mov, dirs)) in shells.iter().zip(&mut cur_dirs) {
            if record_each {
                let prev = dirs[..*n_mov].to_vec();
                moved = moved.max(cvt.step_dirs_pinned(dirs, *n_mov, *radius));
                for (range, (&from, &to)) in shell.iter().zip(prev.iter().zip(dirs.iter())) {
                    rotate_lipid_from_to(positions, range.clone(), center, from, to);
                }
            } else {
                moved = moved.max(cvt.step_dirs_pinned(dirs, *n_mov, *radius));
            }
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
        if moved < converged {
            break;
        }
    }
    if !record_each {
        for ((shell, _), (init, (_, cur))) in shells.iter().zip(init_dirs.iter().zip(&cur_dirs)) {
            for (range, (&from, &to)) in shell.iter().zip(init.iter().zip(cur)) {
                rotate_lipid_from_to(positions, range.clone(), center, from, to);
            }
        }
    }
    relax_on_sphere(
        positions, radii, &all, obstacles, center, iterations, tolerance, observer,
    )
}

/// De-clash `lipids` over their shells so no bead overlaps another lipid's bead
/// or any fixed `obstacles`. Each lipid slides rigidly over its shell about
/// `center` (its bead radii from the centre are preserved).
#[allow(clippy::too_many_arguments)]
fn relax_on_sphere(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    lipids: &[Range<usize>],
    obstacles: &[Obstacle],
    center: [f64; 3],
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let max_radius = radii
        .iter()
        .copied()
        .chain(obstacles.iter().map(|o| o.radius))
        .fold(0.0_f64, f64::max);
    if max_radius <= 0.0 || lipids.is_empty() {
        return 0.0;
    }
    let n_beads = positions.len();
    let n_total = n_beads + obstacles.len();
    let owner = bead_owners(n_beads, lipids);
    let cell = 2.0 * max_radius;
    let inv_cell = 1.0 / cell;
    let max_step = MAX_STEP_FRACTION * max_radius;

    // A flat CSR 3D cell list over the padded bounding box of the lipids + obstacles.
    let pad = 2.0 * cell;
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for p in positions.iter() {
        lo = [lo[0].min(p[0]), lo[1].min(p[1]), lo[2].min(p[2])];
        hi = [hi[0].max(p[0]), hi[1].max(p[1]), hi[2].max(p[2])];
    }
    for o in obstacles {
        lo = [
            lo[0].min(o.pos[0]),
            lo[1].min(o.pos[1]),
            lo[2].min(o.pos[2]),
        ];
        hi = [
            hi[0].max(o.pos[0]),
            hi[1].max(o.pos[1]),
            hi[2].max(o.pos[2]),
        ];
    }
    lo = [lo[0] - pad, lo[1] - pad, lo[2] - pad];
    hi = [hi[0] + pad, hi[1] + pad, hi[2] + pad];
    let dim = |a: usize| (((hi[a] - lo[a]) * inv_cell).ceil() as usize).max(1);
    let (nx, ny, nz) = (dim(0), dim(1), dim(2));
    let ncells = nx * ny * nz;

    let obs_cells: Vec<([f64; 3], [u16; 3], f64)> = obstacles
        .iter()
        .map(|o| {
            let ix = (((o.pos[0] - lo[0]) * inv_cell) as isize).clamp(0, nx as isize - 1) as u16;
            let iy = (((o.pos[1] - lo[1]) * inv_cell) as isize).clamp(0, ny as isize - 1) as u16;
            let iz = (((o.pos[2] - lo[2]) * inv_cell) as isize).clamp(0, nz as isize - 1) as u16;
            (o.pos, [ix, iy, iz], o.radius)
        })
        .collect();

    let mut bead_cell = vec![[0u16; 3]; n_beads];
    let mut counts = vec![0u32; ncells];
    let mut offsets = vec![0u32; ncells + 1];
    let mut sorted = vec![PackBead::default(); n_total];
    let mut force = vec![[0.0_f64; 3]; lipids.len()];

    let mut max_overlap = 0.0_f64;
    let mut tracker = StallTracker::new();
    for _ in 0..iterations {
        counts.fill(0);
        for (i, &p) in positions.iter().enumerate() {
            let ix = (((p[0] - lo[0]) * inv_cell) as isize).clamp(0, nx as isize - 1) as usize;
            let iy = (((p[1] - lo[1]) * inv_cell) as isize).clamp(0, ny as isize - 1) as usize;
            let iz = (((p[2] - lo[2]) * inv_cell) as isize).clamp(0, nz as isize - 1) as usize;
            bead_cell[i] = [ix as u16, iy as u16, iz as u16];
            counts[(iz * ny + iy) * nx + ix] += 1;
        }
        for &(_, [ix, iy, iz], _) in &obs_cells {
            counts[(iz as usize * ny + iy as usize) * nx + ix as usize] += 1;
        }
        prefix_sum(&mut counts, &mut offsets);
        for i in 0..n_beads {
            let [ix, iy, iz] = bead_cell[i];
            let c = (iz as usize * ny + iy as usize) * nx + ix as usize;
            let slot = counts[c] as usize;
            counts[c] += 1;
            sorted[slot] = PackBead {
                pos: positions[i],
                r: radii[i],
                lipid: owner[i],
            };
        }
        for &(pos, [ix, iy, iz], r) in &obs_cells {
            let c = (iz as usize * ny + iy as usize) * nx + ix as usize;
            let slot = counts[c] as usize;
            counts[c] += 1;
            sorted[slot] = PackBead {
                pos,
                r,
                lipid: u32::MAX,
            };
        }

        let pos_ref = &*positions;
        let bead_cell_ref = &bead_cell;
        let offsets_ref = &offsets;
        let sorted_ref = &sorted;

        max_overlap = lipids
            .par_iter()
            .enumerate()
            .zip(force.par_iter_mut())
            .with_min_len(PAR_MIN_LEN)
            .map(|((li, range), f_out)| {
                let li_u32 = li as u32;
                let mut fx = 0.0_f64;
                let mut fy = 0.0_f64;
                let mut fz = 0.0_f64;
                let mut local_max = 0.0_f64;
                for i in range.clone() {
                    let pi = pos_ref[i];
                    let ri = radii[i];
                    let [cx, cy, cz] = bead_cell_ref[i];
                    let (cx, cy, cz) = (cx as usize, cy as usize, cz as usize);
                    let x0 = cx.saturating_sub(1);
                    let x1 = (cx + 1).min(nx - 1);
                    let y0 = cy.saturating_sub(1);
                    let y1 = (cy + 1).min(ny - 1);
                    let z0 = cz.saturating_sub(1);
                    let z1 = (cz + 1).min(nz - 1);

                    for z in z0..=z1 {
                        let z_base = z * ny;
                        for y in y0..=y1 {
                            let row = (z_base + y) * nx;
                            // x is the contiguous axis: scan x0..=x1 as one slice.
                            let start = offsets_ref[row + x0] as usize;
                            let end = offsets_ref[row + x1 + 1] as usize;
                            for bj in &sorted_ref[start..end] {
                                if bj.lipid == li_u32 {
                                    continue;
                                }
                                let dx = pi[0] - bj.pos[0];
                                let dy = pi[1] - bj.pos[1];
                                let dz = pi[2] - bj.pos[2];
                                let d2 = dx * dx + dy * dy + dz * dz;
                                let sigma = ri + bj.r;
                                if d2 < sigma * sigma && d2 > 1e-12 {
                                    let dist = d2.sqrt();
                                    let overlap = sigma - dist;
                                    local_max = local_max.max(overlap);
                                    let scale = overlap / dist;
                                    fx += dx * scale;
                                    fy += dy * scale;
                                    fz += dz * scale;
                                }
                            }
                        }
                    }
                }
                *f_out = [fx, fy, fz];
                local_max
            })
            .reduce(|| 0.0_f64, f64::max);

        if tracker.should_stop(max_overlap, tolerance) {
            break;
        }

        // Slide each overlapping lipid tangentially: rotate it about the sphere
        // centre by the tangential part of its net push, keeping every bead on its shell.
        for (range, &f) in lipids.iter().zip(&force) {
            if f[0] == 0.0 && f[1] == 0.0 && f[2] == 0.0 {
                continue;
            }
            let n = range.len() as f64;
            let mut c = [0.0; 3];
            for i in range.clone() {
                c[0] += positions[i][0];
                c[1] += positions[i][1];
                c[2] += positions[i][2];
            }
            let dir = [
                c[0] / n - center[0],
                c[1] / n - center[1],
                c[2] / n - center[2],
            ];
            let r_center = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            if r_center < 1e-9 {
                continue;
            }
            let inv_r = 1.0 / r_center;
            let d = [dir[0] * inv_r, dir[1] * inv_r, dir[2] * inv_r];
            let radial = f[0] * d[0] + f[1] * d[1] + f[2] * d[2];
            let ft = [
                f[0] - radial * d[0],
                f[1] - radial * d[1],
                f[2] - radial * d[2],
            ];
            let ftmag = (ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2]).sqrt();
            if ftmag < 1e-12 {
                continue;
            }
            let arc = (STEP_FRACTION * ftmag).min(max_step);
            let angle = arc * inv_r;
            let mut axis = [
                d[1] * ft[2] - d[2] * ft[1],
                d[2] * ft[0] - d[0] * ft[2],
                d[0] * ft[1] - d[1] * ft[0],
            ];
            let am = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
            if am < 1e-12 {
                continue;
            }
            let inv_am = 1.0 / am;
            axis = [axis[0] * inv_am, axis[1] * inv_am, axis[2] * inv_am];
            let (sin_a, cos_a) = angle.sin_cos();
            for i in range.clone() {
                let v = [
                    positions[i][0] - center[0],
                    positions[i][1] - center[1],
                    positions[i][2] - center[2],
                ];
                let rv = rotate_with_sin_cos(v, axis, sin_a, cos_a);
                positions[i] = [center[0] + rv[0], center[1] + rv[1], center[2] + rv[2]];
            }
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
    }
    max_overlap
}

// ============================ Cylindrical ============================

/// Extract per-lipid unwrapped cylindrical coordinates `(R·θ, z) ∈ [0, 2πR) × [0, length)`.
fn extract_cylinder_sites(
    positions: &[[f64; 3]],
    shell: &[Range<usize>],
    radius: f64,
    length: f64,
) -> Vec<[f64; 2]> {
    let circumference = TAU * radius;
    shell
        .iter()
        .map(|range| {
            let n = range.len() as f64;
            let (mut cx, mut cy) = (0.0, 0.0);
            for k in range.clone() {
                cx += positions[k][0];
                cy += positions[k][1];
            }
            let theta = (cy / n).atan2(cx / n).rem_euclid(TAU);
            let s = wrap(theta * radius, circumference);
            let z = wrap(positions[range.start][2], length);
            [s, z]
        })
        .collect()
}

/// Rotate every bead of `range` around the `z`-axis by `d_theta` and shift `z`
/// by `dz` modulo `length`, preserving every bead's exact cylindrical radius.
#[inline]
fn move_lipid_on_cylinder(
    positions: &mut [[f64; 3]],
    range: Range<usize>,
    d_theta: f64,
    dz: f64,
    length: f64,
) {
    let (sin_a, cos_a) = d_theta.sin_cos();
    for k in range {
        let [x, y, z] = positions[k];
        positions[k] = [
            x * cos_a - y * sin_a,
            x * sin_a + y * cos_a,
            wrap(z + dz, length),
        ];
    }
}

/// Centroidal-Voronoi relax a cylindrical tubule to a uniform areal density, then
/// de-clash.
///
/// The tubule axis is `z` through `(0, 0)` and periodic in `[0, length)`.
/// `shells` pairs each shell's lipid ranges with its head radius `r` (nm); each
/// shell is relaxed as an independent periodic planar CVT on the isometric
/// unrolling `[0, 2πr) × [0, length)`. `positions` holds every lipid bead
/// (mutated in place); `radii` is the de-clash radius per bead. Returns the
/// largest remaining bead overlap (nm).
pub fn equalize_on_cylinder(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    shells: &[(&[Range<usize>], f64)],
    length: f64,
    iterations: usize,
    tolerance: f64,
    observer: Option<StepObserver>,
) -> f64 {
    equalize_on_cylinder_with_obstacles(
        positions,
        radii,
        shells,
        &[],
        length,
        iterations,
        tolerance,
        observer,
    )
}

/// Like [`equalize_on_cylinder`], but treats `obstacles` as fixed spheres that
/// clip nearby cylindrical Voronoi cells and repel overlapping lipid beads.
#[allow(clippy::too_many_arguments)]
pub fn equalize_on_cylinder_with_obstacles(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    shells: &[(&[Range<usize>], f64)],
    obstacles: &[Obstacle],
    length: f64,
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let all: Vec<Range<usize>> = shells.iter().flat_map(|(s, _)| s.iter().cloned()).collect();
    if all.is_empty() || length <= 0.0 {
        return 0.0;
    }
    let converged = shells
        .iter()
        .filter(|(s, r)| !s.is_empty() && *r > 0.0)
        .map(|(s, r)| CVT_CONVERGED_FRACTION * (TAU * r * length / s.len() as f64).sqrt())
        .fold(f64::MAX, f64::min);
    let mut cvt = PlaneCvtWorkspace::new();

    let max_lipid_r = radii.iter().copied().fold(0.0_f64, f64::max);
    let init_sites: Vec<Vec<[f64; 2]>> = shells
        .iter()
        .map(|(s, r)| extract_cylinder_sites(positions, s, *r, length))
        .collect();
    let mut cur_sites: Vec<(usize, Vec<[f64; 2]>)> = shells
        .iter()
        .zip(&init_sites)
        .map(|((shell, radius), sites)| {
            let n_mov = sites.len();
            let mut combined = sites.clone();
            if !obstacles.is_empty() && !shell.is_empty() && *radius > 0.0 {
                let (mut r_lo, mut r_hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for range in *shell {
                    for k in range.clone() {
                        let r = positions[k][0].hypot(positions[k][1]);
                        r_lo = r_lo.min(r);
                        r_hi = r_hi.max(r);
                    }
                }
                let circ = TAU * *radius;
                let raw: Vec<[f64; 2]> = obstacles
                    .iter()
                    .filter_map(|o| {
                        let r = o.pos[0].hypot(o.pos[1]);
                        let reach = o.radius + max_lipid_r;
                        if r > 1e-6 && r + reach >= r_lo && r - reach <= r_hi {
                            let theta = o.pos[1].atan2(o.pos[0]).rem_euclid(TAU);
                            Some([wrap(theta * *radius, circ), wrap(o.pos[2], length)])
                        } else {
                            None
                        }
                    })
                    .collect();
                combined.extend(dedup_plane_sites(raw, [circ, length], PINNED_SITE_SPACING));
            }
            (n_mov, combined)
        })
        .collect();
    let record_each = observer.is_some();

    for _ in 0..CVT_ROUNDS {
        let mut moved = 0.0_f64;
        for ((shell, radius), (n_mov, sites)) in shells.iter().zip(&mut cur_sites) {
            if shell.is_empty() || *radius <= 0.0 {
                continue;
            }
            let box_sz = [TAU * *radius, length];
            moved = moved.max(cvt.step_sites_pinned(sites, *n_mov, box_sz));
            if record_each {
                let inv_r = 1.0 / *radius;
                for (range, &[ds, dz]) in shell.iter().zip(cvt.last_displacements()) {
                    move_lipid_on_cylinder(positions, range.clone(), ds * inv_r, dz, length);
                }
            }
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
        if moved < converged {
            break;
        }
    }
    if !record_each {
        for ((shell, radius), (init, (_, cur))) in
            shells.iter().zip(init_sites.iter().zip(&cur_sites))
        {
            if shell.is_empty() || *radius <= 0.0 {
                continue;
            }
            let inv_r = 1.0 / *radius;
            for (range, (&s0, &s1)) in shell.iter().zip(init.iter().zip(cur)) {
                let d_theta = (s1[0] - s0[0]) * inv_r;
                let dz = s1[1] - s0[1];
                move_lipid_on_cylinder(positions, range.clone(), d_theta, dz, length);
            }
        }
    }
    relax_on_cylinder(
        positions, radii, &all, obstacles, length, iterations, tolerance, observer,
    )
}

/// De-clash `lipids` on concentric cylindrical shells about the `z`-axis,
/// periodic in `z ∈ [0, length)`. Each lipid slides rigidly in `(θ, z)` while
/// preserving every bead's exact radial distance from the axis.
#[allow(clippy::too_many_arguments)]
fn relax_on_cylinder(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    lipids: &[Range<usize>],
    obstacles: &[Obstacle],
    length: f64,
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let max_radius = radii
        .iter()
        .copied()
        .chain(obstacles.iter().map(|o| o.radius))
        .fold(0.0_f64, f64::max);
    if max_radius <= 0.0 || lipids.is_empty() || length <= 0.0 {
        return 0.0;
    }
    let n_beads = positions.len();
    let n_total = n_beads + obstacles.len();
    let owner = bead_owners(n_beads, lipids);
    let cell = 2.0 * max_radius;
    let inv_cell = 1.0 / cell;
    let max_step = MAX_STEP_FRACTION * max_radius;

    // Non-periodic padded xy square spanning the full cylinder cross-section so
    // rotating lipids never leave the cell grid, combined with periodic z.
    let max_r = positions
        .iter()
        .map(|p| p[0].hypot(p[1]))
        .chain(obstacles.iter().map(|o| o.pos[0].hypot(o.pos[1])))
        .fold(0.0_f64, f64::max);
    let pad = 2.0 * cell;
    let lo_xy = -(max_r + pad);
    let span_xy = 2.0 * (max_r + pad);
    let nxy = ((span_xy * inv_cell).ceil() as usize).max(1);
    let nz = ((length / cell).floor() as usize).max(1);
    let inv_wz = nz as f64 / length;
    let half_z = 0.5 * length;
    let dzs: &[isize] = match nz {
        1 => &[0],
        2 => &[0, 1],
        _ => &[-1, 0, 1],
    };
    let ncells = nxy * nxy * nz;

    let obs_cells: Vec<([f64; 3], [u32; 3], f64)> = obstacles
        .iter()
        .map(|o| {
            let wz = wrap(o.pos[2], length);
            let ix = (((o.pos[0] - lo_xy) * inv_cell) as isize).clamp(0, nxy as isize - 1) as u32;
            let iy = (((o.pos[1] - lo_xy) * inv_cell) as isize).clamp(0, nxy as isize - 1) as u32;
            let iz = ((wz * inv_wz) as usize).min(nz - 1) as u32;
            ([o.pos[0], o.pos[1], wz], [ix, iy, iz], o.radius)
        })
        .collect();

    let mut wrapped = vec![[0.0_f64; 3]; n_beads];
    let mut bead_cell = vec![[0u32; 3]; n_beads];
    let mut counts = vec![0u32; ncells];
    let mut offsets = vec![0u32; ncells + 1];
    let mut sorted = vec![PackBead::default(); n_total];
    let mut force = vec![[0.0_f64; 3]; lipids.len()];

    let mut max_overlap = 0.0_f64;
    let mut tracker = StallTracker::new();
    for _ in 0..iterations {
        counts.fill(0);
        for (i, &p) in positions.iter().enumerate() {
            let wz = wrap(p[2], length);
            let ix = (((p[0] - lo_xy) * inv_cell) as isize).clamp(0, nxy as isize - 1) as usize;
            let iy = (((p[1] - lo_xy) * inv_cell) as isize).clamp(0, nxy as isize - 1) as usize;
            let iz = ((wz * inv_wz) as usize).min(nz - 1);
            wrapped[i] = [p[0], p[1], wz];
            bead_cell[i] = [ix as u32, iy as u32, iz as u32];
            counts[(iz * nxy + iy) * nxy + ix] += 1;
        }
        for &(_, [ix, iy, iz], _) in &obs_cells {
            counts[(iz as usize * nxy + iy as usize) * nxy + ix as usize] += 1;
        }
        prefix_sum(&mut counts, &mut offsets);
        for i in 0..n_beads {
            let [ix, iy, iz] = bead_cell[i];
            let c = (iz as usize * nxy + iy as usize) * nxy + ix as usize;
            let slot = counts[c] as usize;
            counts[c] += 1;
            sorted[slot] = PackBead {
                pos: wrapped[i],
                r: radii[i],
                lipid: owner[i],
            };
        }
        for &(pos, [ix, iy, iz], r) in &obs_cells {
            let c = (iz as usize * nxy + iy as usize) * nxy + ix as usize;
            let slot = counts[c] as usize;
            counts[c] += 1;
            sorted[slot] = PackBead {
                pos,
                r,
                lipid: u32::MAX,
            };
        }

        let wrapped_ref = &wrapped;
        let bead_cell_ref = &bead_cell;
        let offsets_ref = &offsets;
        let sorted_ref = &sorted;

        max_overlap = lipids
            .par_iter()
            .enumerate()
            .zip(force.par_iter_mut())
            .with_min_len(PAR_MIN_LEN)
            .map(|((li, range), f_out)| {
                let li_u32 = li as u32;
                let mut fx = 0.0_f64;
                let mut fy = 0.0_f64;
                let mut fz = 0.0_f64;
                let mut local_max = 0.0_f64;
                for i in range.clone() {
                    let pi = wrapped_ref[i];
                    let ri = radii[i];
                    let [cx, cy, cz] = bead_cell_ref[i];
                    let x0 = (cx as usize).saturating_sub(1);
                    let x1 = (cx as usize + 1).min(nxy - 1);
                    let y0 = (cy as usize).saturating_sub(1);
                    let y1 = (cy as usize + 1).min(nxy - 1);

                    for &dzc in dzs {
                        let z_base = wrap_cell(cz, dzc, nz) * nxy;
                        for y in y0..=y1 {
                            let row = (z_base + y) * nxy;
                            let start = offsets_ref[row + x0] as usize;
                            let end = offsets_ref[row + x1 + 1] as usize;
                            for bj in &sorted_ref[start..end] {
                                if bj.lipid == li_u32 {
                                    continue;
                                }
                                let dx = pi[0] - bj.pos[0];
                                let dy = pi[1] - bj.pos[1];
                                let mut dz = pi[2] - bj.pos[2];
                                if dz > half_z {
                                    dz -= length;
                                } else if dz < -half_z {
                                    dz += length;
                                }
                                let d2 = dx * dx + dy * dy + dz * dz;
                                let sigma = ri + bj.r;
                                if d2 < sigma * sigma && d2 > 1e-12 {
                                    let dist = d2.sqrt();
                                    let overlap = sigma - dist;
                                    local_max = local_max.max(overlap);
                                    let scale = overlap / dist;
                                    fx += dx * scale;
                                    fy += dy * scale;
                                    fz += dz * scale;
                                }
                            }
                        }
                    }
                }
                *f_out = [fx, fy, fz];
                local_max
            })
            .reduce(|| 0.0_f64, f64::max);

        if tracker.should_stop(max_overlap, tolerance) {
            break;
        }

        // Slide each overlapping lipid tangentially in θ and axially in z,
        // preserving every bead's exact radial distance from the cylinder axis.
        for (range, &[fx, fy, fz]) in lipids.iter().zip(&force) {
            if fx == 0.0 && fy == 0.0 && fz == 0.0 {
                continue;
            }
            let n = range.len() as f64;
            let (mut cx, mut cy) = (0.0, 0.0);
            for i in range.clone() {
                cx += positions[i][0];
                cy += positions[i][1];
            }
            cx /= n;
            cy /= n;
            let r_c = cx.hypot(cy);
            if r_c < 1e-9 {
                continue;
            }
            let inv_r = 1.0 / r_c;
            let tx = -cy * inv_r;
            let ty = cx * inv_r;
            let f_theta = fx * tx + fy * ty;
            let f_mag = f_theta.hypot(fz);
            if f_mag < 1e-12 {
                continue;
            }
            let step = (STEP_FRACTION * f_mag).min(max_step);
            let ds = step * (f_theta / f_mag);
            let dz = step * (fz / f_mag);
            move_lipid_on_cylinder(positions, range.clone(), ds * inv_r, dz, length);
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
    }
    max_overlap
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relax_separates_two_overlapping_lipids() {
        // Two single-bead lipids overlapping (centres 0.5 nm apart, sigma 1.0).
        let mut pos = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let radii = [0.5, 0.5];
        let lipids = [0..1, 1..2];
        let max = relax(
            &mut pos,
            &radii,
            &lipids,
            &[],
            [20.0, 20.0],
            500,
            1e-3,
            None,
        );
        let d = ((pos[0][0] - pos[1][0]).powi(2) + (pos[0][1] - pos[1][1]).powi(2)).sqrt();
        assert!(max < 1e-3, "residual overlap {max}");
        assert!(d >= 1.0 - 1e-2, "beads still overlap: d = {d}");
    }

    #[test]
    fn equalize_plane_with_obstacles_keeps_lipids_outside_obstacle() {
        let box_xy = [16.0, 16.0];
        let obs = Obstacle {
            pos: [8.0, 8.0, 0.0],
            radius: 2.0,
        };
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        for i in 0..10 {
            for j in 0..10 {
                let x = 0.8 + i as f64 * 1.5;
                let y = 0.8 + j as f64 * 1.5;
                if (x - 8.0).hypot(y - 8.0) < 2.6 {
                    continue;
                }
                let k = positions.len();
                positions.push([x, y, 0.0]);
                lipids.push(k..k + 1);
            }
        }
        let radii = vec![0.4; positions.len()];
        let leaflet: Vec<_> = lipids.clone();
        let overlap = equalize_plane_with_obstacles(
            &mut positions,
            &radii,
            &[&leaflet],
            &[obs],
            box_xy,
            1000,
            0.01,
            None,
        );
        assert!(overlap < 0.05, "residual obstacle overlap {overlap}");
        for p in &positions {
            let dx = (p[0].rem_euclid(16.0) - 8.0)
                .abs()
                .min(16.0 - (p[0].rem_euclid(16.0) - 8.0).abs());
            let dy = (p[1].rem_euclid(16.0) - 8.0)
                .abs()
                .min(16.0 - (p[1].rem_euclid(16.0) - 8.0).abs());
            assert!(
                dx.hypot(dy) >= 2.35,
                "lipid penetrated obstacle: d = {}",
                dx.hypot(dy)
            );
        }
    }

    #[test]
    fn equalize_plane_leaves_a_clash_free_uniform_leaflet() {
        // 256 single-bead lipids placed in one corner of a 32×32 periodic box: the
        // CVT must flow them out to fill the box, the de-clash must remove overlaps.
        let box_xy = [32.0, 32.0];
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        for i in 0..16 {
            for j in 0..16 {
                let k = positions.len();
                positions.push([0.2 * i as f64, 0.2 * j as f64, 0.0]);
                lipids.push(k..k + 1);
            }
        }
        let radii = vec![0.5; positions.len()];
        let leaflet: Vec<_> = lipids.clone();
        let overlap = equalize_plane(
            &mut positions,
            &radii,
            &[&leaflet],
            box_xy,
            1000,
            0.01,
            None,
        );
        assert!(overlap < 0.05, "not clash-free: max overlap {overlap}");
        // The lipids now span the box rather than a corner.
        let span = |axis: usize| {
            let mut lo = f64::MAX;
            let mut hi = f64::MIN;
            for p in &positions {
                let v = p[axis].rem_euclid(box_xy[axis]);
                lo = lo.min(v);
                hi = hi.max(v);
            }
            hi - lo
        };
        assert!(span(0) > 24.0 && span(1) > 24.0, "leaflet did not spread");
    }

    #[test]
    fn equalize_on_sphere_is_clash_free_and_stays_on_the_shell() {
        // Single-bead lipids clustered on a cap of a radius-10 sphere; the spherical
        // CVT spreads them, the de-clash removes overlaps, all stay on the shell.
        let center = [0.0, 0.0, 0.0];
        let r = 10.0;
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        for i in 0..12 {
            for j in 0..12 {
                let a = (i as f64 - 5.5) * 0.06;
                let b = (j as f64 - 5.5) * 0.06;
                let m = (a * a + b * b + 1.0_f64).sqrt();
                let k = positions.len();
                positions.push([a / m * r, b / m * r, 1.0 / m * r]);
                lipids.push(k..k + 1);
            }
        }
        let radii = vec![0.5; positions.len()];
        let shell: Vec<_> = lipids.clone();
        let overlap = equalize_on_sphere(
            &mut positions,
            &radii,
            &[(&shell, r)],
            center,
            1000,
            0.01,
            None,
        );
        assert!(overlap < 0.1, "not clash-free: max overlap {overlap}");
        for p in &positions {
            let radius = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(
                (radius - r).abs() < 1e-6,
                "bead left the shell: r = {radius}"
            );
        }
    }

    #[test]
    fn equalize_on_cylinder_is_clash_free_and_preserves_radius() {
        // Single-bead lipids clustered on one patch of a cylinder of radius 8,
        // length 20; cylindrical CVT must spread them around θ and z while
        // de-clash removes overlaps and preserves every bead's exact radius.
        let r = 8.0;
        let length = 20.0;
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        for i in 0..12 {
            for j in 0..12 {
                let theta = i as f64 * 0.04;
                let z = 1.0 + j as f64 * 0.3;
                let k = positions.len();
                positions.push([r * theta.cos(), r * theta.sin(), z]);
                lipids.push(k..k + 1);
            }
        }
        let radii = vec![0.5; positions.len()];
        let shell: Vec<_> = lipids.clone();
        let overlap = equalize_on_cylinder(
            &mut positions,
            &radii,
            &[(&shell, r)],
            length,
            1000,
            0.01,
            None,
        );
        assert!(overlap < 0.05, "not clash-free: max overlap {overlap}");
        let mut z_lo = f64::MAX;
        let mut z_hi = f64::MIN;
        let mut x_min = f64::MAX;
        for p in &positions {
            let radius = p[0].hypot(p[1]);
            assert!(
                (radius - r).abs() < 1e-6,
                "bead left cylinder shell: r = {radius}"
            );
            z_lo = z_lo.min(p[2]);
            z_hi = z_hi.max(p[2]);
            x_min = x_min.min(p[0]);
        }
        assert!(z_hi - z_lo > 14.0, "lipids did not spread axially");
        assert!(x_min < -4.0, "lipids did not spread circumferentially");
    }
}
