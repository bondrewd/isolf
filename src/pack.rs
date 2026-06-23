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
//!
//! Beads are treated as soft spheres: an overlapping pair pushes apart by a fraction
//! of the overlap each step (a damped Jacobi relaxation), capped for stability, over
//! a cell list. A whole lipid moves rigidly (laterally in the plane, with z fixed so
//! it stays in its leaflet, or tangentially over the sphere at fixed radius), so the
//! relaxation never distorts a lipid or moves it out of its layer. This is the clean
//! starting point; an NPT/NVT run then finds the equilibrium density.

use std::ops::Range;

use crate::voronoi::{cvt_plane_step, cvt_sphere_step, rotate_about_axis};

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

// ============================ Planar ============================

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
    mut observer: Option<StepObserver>,
) -> f64 {
    let all: Vec<Range<usize>> = leaflets.iter().flat_map(|l| l.iter().cloned()).collect();
    let total: usize = leaflets.iter().map(|l| l.len()).sum();
    if total == 0 || box_xy[0] <= 0.0 || box_xy[1] <= 0.0 {
        return 0.0;
    }
    let converged = CVT_CONVERGED_FRACTION * (box_xy[0] * box_xy[1] / total as f64).sqrt();
    for _ in 0..CVT_ROUNDS {
        let mut moved = 0.0_f64;
        for leaflet in leaflets {
            moved = moved.max(cvt_plane_step(positions, leaflet, box_xy));
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
        if moved < converged {
            break;
        }
    }
    relax(
        positions, radii, &all, box_xy, iterations, tolerance, observer,
    )
}

/// De-clash `lipids` laterally so no bead overlaps another lipid's bead, in a
/// periodic box `box_xy` (nm). Each lipid translates rigidly in x/y; z is fixed.
/// Stops when the largest overlap drops below `tolerance` (nm) or after
/// `iterations`. Returns the largest remaining overlap.
fn relax(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    lipids: &[Range<usize>],
    box_xy: [f64; 2],
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let mut owner = vec![usize::MAX; positions.len()];
    for (lipid, range) in lipids.iter().enumerate() {
        for i in range.clone() {
            owner[i] = lipid;
        }
    }
    let max_radius = radii.iter().copied().fold(0.0_f64, f64::max);
    if max_radius <= 0.0 || lipids.is_empty() {
        return 0.0;
    }
    // Cells at least one bead diameter wide, so any overlapping pair is in adjacent
    // cells. The grid covers the periodic x/y box.
    let cell = 2.0 * max_radius;
    let nx = ((box_xy[0] / cell).floor() as usize).max(1);
    let ny = ((box_xy[1] / cell).floor() as usize).max(1);
    let (wx, wy) = (box_xy[0] / nx as f64, box_xy[1] / ny as f64);
    let cell_of = |x: f64, y: f64| {
        let gx = ((x.rem_euclid(box_xy[0])) / wx).floor() as usize % nx;
        let gy = ((y.rem_euclid(box_xy[1])) / wy).floor() as usize % ny;
        (gx, gy)
    };
    // Neighbour-cell offsets per axis. With one or two cells a ±1 sweep would wrap
    // back onto a cell already visited and double-count its beads, so list each
    // reachable offset exactly once.
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

    // Scratch reused across iterations so the hot loop never reallocates. Packing
    // each bead's position, radius, and owner into one tuple keeps a neighbour read
    // on a single cache line, which outweighs the per-step copy on dense systems
    // (measured: indexing the three slices directly is ~5% slower on a vesicle).
    let mut beads: Vec<([f64; 3], f64, usize)> = Vec::with_capacity(positions.len());
    let mut grid: Vec<Vec<usize>> = vec![Vec::new(); nx * ny];
    let mut force = vec![[0.0_f64; 2]; lipids.len()];
    let mut max_overlap = 0.0_f64;
    let mut best = f64::MAX;
    let mut stalled = 0usize;
    for _ in 0..iterations {
        beads.clear();
        beads.extend(
            positions
                .iter()
                .zip(radii)
                .enumerate()
                .map(|(i, (&p, &r))| (p, r, owner[i])),
        );

        for bucket in grid.iter_mut() {
            bucket.clear();
        }
        for (i, &(p, _, _)) in beads.iter().enumerate() {
            let (gx, gy) = cell_of(p[0], p[1]);
            grid[gy * nx + gx].push(i);
        }

        for f in force.iter_mut() {
            *f = [0.0; 2];
        }
        max_overlap = 0.0;
        for (i, &(pi, ri, li)) in beads.iter().enumerate() {
            let (gx, gy) = cell_of(pi[0], pi[1]);
            for &dyc in dys {
                for &dxc in dxs {
                    let cx = (gx as isize + dxc).rem_euclid(nx as isize) as usize;
                    let cy = (gy as isize + dyc).rem_euclid(ny as isize) as usize;
                    for &j in &grid[cy * nx + cx] {
                        if j == i {
                            continue;
                        }
                        let (pj, rj, lj) = beads[j];
                        if li == lj {
                            continue; // same lipid: rigid, no internal force
                        }
                        let mut dx = pi[0] - pj[0];
                        let mut dy = pi[1] - pj[1];
                        dx -= box_xy[0] * (dx / box_xy[0]).round();
                        dy -= box_xy[1] * (dy / box_xy[1]).round();
                        let dz = pi[2] - pj[2];
                        let d2 = dx * dx + dy * dy + dz * dz;
                        let sigma = ri + rj;
                        if d2 < sigma * sigma && d2 > 1e-12 {
                            let d = d2.sqrt();
                            let overlap = sigma - d;
                            max_overlap = max_overlap.max(overlap);
                            force[li][0] += overlap * dx / d;
                            force[li][1] += overlap * dy / d;
                        }
                    }
                }
            }
        }

        if max_overlap < tolerance {
            break;
        }
        if best - max_overlap < MIN_IMPROVEMENT {
            stalled += 1;
            if stalled >= STALL_LIMIT {
                break;
            }
        } else {
            stalled = 0;
        }
        best = best.min(max_overlap);

        for (lipid, range) in lipids.iter().enumerate() {
            let mut sx = STEP_FRACTION * force[lipid][0];
            let mut sy = STEP_FRACTION * force[lipid][1];
            let mag = (sx * sx + sy * sy).sqrt();
            if mag > max_step {
                sx *= max_step / mag;
                sy *= max_step / mag;
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
    for _ in 0..CVT_ROUNDS {
        let mut moved = 0.0_f64;
        for (shell, radius) in shells {
            moved = moved.max(cvt_sphere_step(positions, shell, center, *radius));
        }
        if let Some(observe) = observer.as_deref_mut() {
            observe(positions);
        }
        if moved < converged {
            break;
        }
    }
    relax_on_sphere(
        positions, radii, &all, center, iterations, tolerance, observer,
    )
}

/// De-clash `lipids` over their shells so no bead overlaps another lipid's bead.
/// Each lipid slides rigidly over its shell about `center` (its bead radii from the
/// centre are preserved). Stops when the largest overlap drops below `tolerance`
/// (nm) or after `iterations`. Returns the largest remaining overlap.
fn relax_on_sphere(
    positions: &mut [[f64; 3]],
    radii: &[f64],
    lipids: &[Range<usize>],
    center: [f64; 3],
    iterations: usize,
    tolerance: f64,
    mut observer: Option<StepObserver>,
) -> f64 {
    let mut owner = vec![usize::MAX; positions.len()];
    for (lipid, range) in lipids.iter().enumerate() {
        for i in range.clone() {
            owner[i] = lipid;
        }
    }
    let max_radius = radii.iter().copied().fold(0.0_f64, f64::max);
    if max_radius <= 0.0 || lipids.is_empty() {
        return 0.0;
    }
    let cell = 2.0 * max_radius;
    let max_step = MAX_STEP_FRACTION * max_radius;

    // A dense 3D cell list over the bounding box of the lipids. Cells one bead
    // diameter wide put any overlap in an adjacent cell; the box is padded so a
    // sliding lipid stays inside, and indices clamp.
    let pad = 2.0 * cell;
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for p in positions.iter() {
        lo = [lo[0].min(p[0]), lo[1].min(p[1]), lo[2].min(p[2])];
        hi = [hi[0].max(p[0]), hi[1].max(p[1]), hi[2].max(p[2])];
    }
    lo = [lo[0] - pad, lo[1] - pad, lo[2] - pad];
    hi = [hi[0] + pad, hi[1] + pad, hi[2] + pad];
    let dim = |a: usize| (((hi[a] - lo[a]) / cell).ceil() as usize).max(1);
    let (nx, ny, nz) = (dim(0), dim(1), dim(2));
    let coords = |p: [f64; 3]| {
        let ax = |v: f64, l: f64, n: usize| {
            (((v - l) / cell) as isize).clamp(0, n as isize - 1) as usize
        };
        (
            ax(p[0], lo[0], nx),
            ax(p[1], lo[1], ny),
            ax(p[2], lo[2], nz),
        )
    };
    let mut grid: Vec<Vec<usize>> = vec![Vec::new(); nx * ny * nz];
    // Beads packed (position, radius, owner) so each neighbour read is one cache
    // line; see the note in `relax`.
    let mut beads: Vec<([f64; 3], f64, usize)> = Vec::with_capacity(positions.len());
    let mut force = vec![[0.0_f64; 3]; lipids.len()];

    let mut max_overlap = 0.0_f64;
    let mut best = f64::MAX;
    let mut stalled = 0usize;
    for _ in 0..iterations {
        beads.clear();
        beads.extend(
            positions
                .iter()
                .zip(radii)
                .enumerate()
                .map(|(i, (&p, &r))| (p, r, owner[i])),
        );

        for bucket in grid.iter_mut() {
            bucket.clear();
        }
        for (i, &(p, _, _)) in beads.iter().enumerate() {
            let (ix, iy, iz) = coords(p);
            grid[ix + nx * (iy + ny * iz)].push(i);
        }

        for f in force.iter_mut() {
            *f = [0.0; 3];
        }
        max_overlap = 0.0;
        for (i, &(pi, ri, li)) in beads.iter().enumerate() {
            let (cx, cy, cz) = coords(pi);
            for dz in -1..=1_isize {
                let z = cz as isize + dz;
                if z < 0 || z >= nz as isize {
                    continue;
                }
                for dy in -1..=1_isize {
                    let y = cy as isize + dy;
                    if y < 0 || y >= ny as isize {
                        continue;
                    }
                    for dx in -1..=1_isize {
                        let x = cx as isize + dx;
                        if x < 0 || x >= nx as isize {
                            continue;
                        }
                        let bucket = &grid[x as usize + nx * (y as usize + ny * z as usize)];
                        for &j in bucket {
                            if j == i {
                                continue;
                            }
                            let (pj, rj, lj) = beads[j];
                            if li == lj {
                                continue; // same lipid: rigid, no internal force
                            }
                            let d = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                            let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                            let sigma = ri + rj;
                            if d2 < sigma * sigma && d2 > 1e-12 {
                                let dist = d2.sqrt();
                                let overlap = sigma - dist;
                                max_overlap = max_overlap.max(overlap);
                                for a in 0..3 {
                                    force[li][a] += overlap * d[a] / dist;
                                }
                            }
                        }
                    }
                }
            }
        }

        if max_overlap < tolerance {
            break;
        }
        if best - max_overlap < MIN_IMPROVEMENT {
            stalled += 1;
            if stalled >= STALL_LIMIT {
                break;
            }
        } else {
            stalled = 0;
        }
        best = best.min(max_overlap);

        // Slide each lipid tangentially: rotate it about the sphere centre by the
        // tangential part of its net push, keeping every bead on its shell.
        for (lipid, range) in lipids.iter().enumerate() {
            let n = range.len() as f64;
            let mut c = [0.0; 3];
            for i in range.clone() {
                for a in 0..3 {
                    c[a] += positions[i][a];
                }
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
            let d = [dir[0] / r_center, dir[1] / r_center, dir[2] / r_center];
            let f = force[lipid];
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
            let angle = arc / r_center;
            let mut axis = [
                d[1] * ft[2] - d[2] * ft[1],
                d[2] * ft[0] - d[0] * ft[2],
                d[0] * ft[1] - d[1] * ft[0],
            ];
            let am = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
            if am < 1e-12 {
                continue;
            }
            for x in &mut axis {
                *x /= am;
            }
            for i in range.clone() {
                let v = [
                    positions[i][0] - center[0],
                    positions[i][1] - center[1],
                    positions[i][2] - center[2],
                ];
                let rv = rotate_about_axis(v, axis, angle);
                positions[i] = [center[0] + rv[0], center[1] + rv[1], center[2] + rv[2]];
            }
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
        let max = relax(&mut pos, &radii, &lipids, [20.0, 20.0], 500, 1e-3, None);
        let d = ((pos[0][0] - pos[1][0]).powi(2) + (pos[0][1] - pos[1][1]).powi(2)).sqrt();
        assert!(max < 1e-3, "residual overlap {max}");
        assert!(d >= 1.0 - 1e-2, "beads still overlap: d = {d}");
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
}
