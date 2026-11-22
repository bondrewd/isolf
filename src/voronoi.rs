//! Centroidal Voronoi tessellation (Lloyd's relaxation) by **per-site clipping**.
//!
//! A site's Voronoi cell is, by definition, the intersection of the half-spaces
//! `{x : x is closer to the site than to neighbour n}` over all neighbours `n`.
//! Rather than build a global Delaunay/convex-hull triangulation, each cell is
//! computed independently by clipping a starting polygon against the perpendicular
//! bisectors of nearby sites (Sutherland–Hodgman), ordered nearest-first with a
//! two-stage adaptive radius search and exact circumradius early termination.
//! This yields the exact Voronoi cell, is robust (no global topology to get
//! wrong), parallelises per site, and needs no external geometry dependency.
//!
//! One Lloyd step moves each site to its cell's centroid; iterating drives the
//! sites toward a uniform (centroidal) distribution. Two surfaces are supported:
//!
//! - **Plane** (a flat bilayer leaflet): a periodic box, bisector half-planes,
//!   polygon area centroid. [`cvt_plane_step`].
//! - **Sphere** (a vesicle shell): bisector planes through the origin (great-circle
//!   half-spaces), spherical-polygon area-weighted centroid. [`cvt_sphere_step`].
//!
//! Each "site" is one lipid, represented by its lateral centroid; the whole lipid
//! is moved rigidly (a translation in the plane, a rotation about the centre on the
//! sphere) so every bead follows and the lipid stays in its leaflet/shell.

use std::ops::Range;

use rayon::prelude::*;

/// Outer neighbour-search radius, in mean lipid spacings: large enough to bound a
/// cell even next to voids in the initial random placement.
const OUTER_SPACINGS: f64 = 4.0;

/// Inner neighbour-search radius, in mean lipid spacings: covers all Voronoi
/// neighbours of >99% of sites once spread; sites whose polygon circumradius
/// could still be clipped by a farther point fall back to [`OUTER_SPACINGS`].
const INNER_SPACINGS: f64 = 2.0;

/// Largest cap/neighbour angular radius (rad), so a small lipid count (huge `beta`)
/// can't ask for a cap near or past a hemisphere.
const MAX_CAP_ALPHA: f64 = 1.0;

/// Precomputed `(cos θ, sin θ)` for the 12-gon starting cap polygon on the sphere
/// (`θ = 2πk/12` for `k = 0..12`). Avoids trig calls per site.
const CAP_CS: [(f64, f64); 12] = [
    (1.0, 0.0),
    (0.866_025_403_784_438_7, 0.5),
    (0.5, 0.866_025_403_784_438_7),
    (0.0, 1.0),
    (-0.5, 0.866_025_403_784_438_7),
    (-0.866_025_403_784_438_7, 0.5),
    (-1.0, 0.0),
    (-0.866_025_403_784_438_7, -0.5),
    (-0.5, -0.866_025_403_784_438_7),
    (0.0, -1.0),
    (0.5, -0.866_025_403_784_438_7),
    (0.866_025_403_784_438_7, -0.5),
];

/// Minimum chunk size for Rayon parallel site sweeps so tiny systems avoid
/// thread-pool dispatch overhead.
const PAR_MIN_LEN: usize = 64;

/// Fast periodic wrap into `[0, len)` for values already near the box.
#[inline]
pub(crate) fn wrap(v: f64, len: f64) -> f64 {
    if v >= 0.0 && v < len {
        v
    } else if v >= len && v < 2.0 * len {
        v - len
    } else if v < 0.0 && v >= -len {
        v + len
    } else {
        v.rem_euclid(len)
    }
}

// ============================ Planar ============================

/// Clip convex polygon `poly` (in site-relative coordinates) into `out`, keeping
/// `{p : p·dir ≤ rhs}` where `dir` points from the site to a neighbour and
/// `rhs = 0.5 * |dir|²`. Returns the squared circumradius `max |p|²` of `out`.
fn clip_halfplane_rel(poly: &[[f64; 2]], dir: [f64; 2], rhs: f64, out: &mut Vec<[f64; 2]>) -> f64 {
    out.clear();
    let n = poly.len();
    if n == 0 {
        return 0.0;
    }
    let mut r_max_sq = 0.0_f64;
    let mut a = poly[n - 1];
    let mut sa = a[0] * dir[0] + a[1] * dir[1] - rhs;
    let mut a_in = sa <= 0.0;
    for &b in poly {
        let sb = b[0] * dir[0] + b[1] * dir[1] - rhs;
        let b_in = sb <= 0.0;
        if a_in != b_in {
            let t = sa / (sa - sb);
            let p = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])];
            r_max_sq = r_max_sq.max(p[0] * p[0] + p[1] * p[1]);
            out.push(p);
        }
        if b_in {
            r_max_sq = r_max_sq.max(b[0] * b[0] + b[1] * b[1]);
            out.push(b);
        }
        a = b;
        sa = sb;
        a_in = b_in;
    }
    r_max_sq
}

/// Area centroid of a simple polygon (the shoelace formula). `None` when the
/// polygon is degenerate (near-zero area).
fn polygon_centroid(poly: &[[f64; 2]]) -> Option<[f64; 2]> {
    let n = poly.len();
    if n < 3 {
        return None;
    }
    let (mut a2, mut cx, mut cy) = (0.0, 0.0, 0.0);
    let mut p = poly[n - 1];
    for &q in poly {
        let cross = p[0] * q[1] - q[0] * p[1];
        a2 += cross;
        cx += (p[0] + q[0]) * cross;
        cy += (p[1] + q[1]) * cross;
        p = q;
    }
    if a2.abs() < 1e-12 {
        return None;
    }
    Some([cx / (3.0 * a2), cy / (3.0 * a2)])
}

/// Clip `poly` against distance-sorted relative neighbours `(d2, dir)`, stopping
/// once `d2 > 4 * r_max_sq`. Updates `r_max_sq` in place; returns `false` if the
/// cell collapsed below 3 vertices.
fn clip_sorted_plane(
    neighbours: &[(f64, [f64; 2])],
    poly: &mut Vec<[f64; 2]>,
    next: &mut Vec<[f64; 2]>,
    r_max_sq: &mut f64,
) -> bool {
    for &(d2, dir) in neighbours {
        if d2 > 4.0 * *r_max_sq * (1.0 + 1e-12) {
            break;
        }
        *r_max_sq = clip_halfplane_rel(poly, dir, 0.5 * d2, next);
        std::mem::swap(poly, next);
        if poly.len() < 3 {
            return false;
        }
    }
    true
}

#[derive(Clone, Copy, Default)]
struct PlaneSite {
    x: f64,
    y: f64,
    id: u32,
}

struct PlaneScratch {
    neighbours: Vec<(f64, [f64; 2])>,
    poly: Vec<[f64; 2]>,
    next: Vec<[f64; 2]>,
}

impl PlaneScratch {
    fn new() -> Self {
        Self {
            neighbours: Vec::with_capacity(32),
            poly: Vec::with_capacity(16),
            next: Vec::with_capacity(16),
        }
    }
}

/// Extract the lateral `(x, y)` centroid of each lipid, wrapped into `[0, box_xy)`.
pub(crate) fn extract_plane_sites(
    positions: &[[f64; 3]],
    lipids: &[Range<usize>],
    box_xy: [f64; 2],
) -> Vec<[f64; 2]> {
    lipids
        .iter()
        .map(|range| {
            let m = range.len() as f64;
            let (mut x, mut y) = (0.0, 0.0);
            for k in range.clone() {
                x += positions[k][0];
                y += positions[k][1];
            }
            [wrap(x / m, box_xy[0]), wrap(y / m, box_xy[1])]
        })
        .collect()
}

/// Reusable scratch buffers for multi-round planar CVT (zero allocations across
/// rounds).
pub(crate) struct PlaneCvtWorkspace {
    cell_coords: Vec<(u32, u32)>,
    counts: Vec<u32>,
    offsets: Vec<u32>,
    sorted: Vec<PlaneSite>,
    disps: Vec<[f64; 2]>,
}

impl PlaneCvtWorkspace {
    pub(crate) fn new() -> Self {
        Self {
            cell_coords: Vec::new(),
            counts: Vec::new(),
            offsets: Vec::new(),
            sorted: Vec::new(),
            disps: Vec::new(),
        }
    }

    /// Advance 2D periodic sites `sites` (already in `[0, box_xy)`) by one exact
    /// Lloyd step in place. Returns the largest in-plane displacement (nm).
    pub(crate) fn step_sites(&mut self, sites: &mut [[f64; 2]], box_xy: [f64; 2]) -> f64 {
        self.step_sites_pinned(sites, sites.len(), box_xy)
    }

    /// Like [`step_sites`](Self::step_sites), but only `sites[..n_movable]`
    /// move; `sites[n_movable..]` stay fixed as obstacle boundary sites that
    /// still clip movable Voronoi cells.
    pub(crate) fn step_sites_pinned(
        &mut self,
        sites: &mut [[f64; 2]],
        n_movable: usize,
        box_xy: [f64; 2],
    ) -> f64 {
        let n = sites.len();
        let n_movable = n_movable.min(n);
        let (bx, by) = (box_xy[0], box_xy[1]);
        if n_movable == 0 || n < 2 || bx <= 0.0 || by <= 0.0 {
            self.disps.resize(n_movable, [0.0; 2]);
            self.disps.fill([0.0; 2]);
            return 0.0;
        }

        let spacing = (bx * by / n as f64).sqrt();
        let r_outer = (OUTER_SPACINGS * spacing).min(0.49 * bx.min(by));
        let r_inner = (INNER_SPACINGS * spacing).min(r_outer);
        let r_inner2 = r_inner * r_inner;
        let r_outer2 = r_outer * r_outer;

        // Cell size >= r_inner so the fast-path 3x3 scan covers r_inner completely.
        let target = r_inner.max(spacing);
        let nx = ((bx / target).floor() as usize).max(1);
        let ny = ((by / target).floor() as usize).max(1);
        let inv_cx = nx as f64 / bx;
        let inv_cy = ny as f64 / by;
        let ncells = nx * ny;

        self.cell_coords.clear();
        self.cell_coords.reserve(n);
        self.counts.resize(ncells, 0);
        self.counts.fill(0);

        for &[wx, wy] in sites.iter() {
            let gx = ((wx * inv_cx) as usize).min(nx - 1);
            let gy = ((wy * inv_cy) as usize).min(ny - 1);
            self.cell_coords.push((gx as u32, gy as u32));
            self.counts[gy * nx + gx] += 1;
        }

        self.offsets.resize(ncells + 1, 0);
        let mut acc = 0u32;
        for c in 0..ncells {
            self.offsets[c] = acc;
            acc += self.counts[c];
            self.counts[c] = self.offsets[c];
        }
        self.offsets[ncells] = acc;

        self.sorted.resize(n, PlaneSite::default());
        for i in 0..n {
            let (gx, gy) = self.cell_coords[i];
            let c = gy as usize * nx + gx as usize;
            let slot = self.counts[c] as usize;
            self.counts[c] += 1;
            self.sorted[slot] = PlaneSite {
                x: sites[i][0],
                y: sites[i][1],
                id: i as u32,
            };
        }

        let rx1 = (r_inner * inv_cx).ceil() as isize;
        let ry1 = (r_inner * inv_cy).ceil() as isize;
        let rx2 = (r_outer * inv_cx).ceil() as isize;
        let ry2 = (r_outer * inv_cy).ceil() as isize;
        let full_x1 = (2 * rx1 + 1) as usize >= nx;
        let full_y1 = (2 * ry1 + 1) as usize >= ny;
        let full_x2 = (2 * rx2 + 1) as usize >= nx;
        let full_y2 = (2 * ry2 + 1) as usize >= ny;
        let half_bx = 0.5 * bx;
        let half_by = 0.5 * by;
        let same_radius = r_inner >= r_outer;

        let sites_ref = &*sites;
        let cell_coords = &self.cell_coords;
        let offsets = &self.offsets;
        let sorted = &self.sorted;

        self.disps.resize(n_movable, [0.0; 2]);
        self.disps
            .par_iter_mut()
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each_init(PlaneScratch::new, |scratch, (i, out)| {
                let [sx, sy] = sites_ref[i];
                let (gx, gy) = (cell_coords[i].0 as isize, cell_coords[i].1 as isize);
                let id_u32 = i as u32;

                // Helper to scan a ring and collect neighbours with min_d2 < d2 <= max_d2.
                let gather = |rx: isize,
                              ry: isize,
                              full_x: bool,
                              full_y: bool,
                              min_d2: f64,
                              max_d2: f64,
                              buf: &mut Vec<(f64, [f64; 2])>| {
                    buf.clear();
                    let nx_i = nx as isize;
                    let ny_i = ny as isize;
                    let y_range = if full_y { 0..ny_i } else { -ry..ry + 1 };
                    let x_range = if full_x { 0..nx_i } else { -rx..rx + 1 };
                    for dy_c in y_range {
                        let cy = if full_y {
                            dy_c as usize
                        } else {
                            let mut c = gy + dy_c;
                            if c < 0 {
                                c += ny_i;
                            } else if c >= ny_i {
                                c -= ny_i;
                            }
                            c as usize
                        };
                        let row = cy * nx;
                        for dx_c in x_range.clone() {
                            let cx = if full_x {
                                dx_c as usize
                            } else {
                                let mut c = gx + dx_c;
                                if c < 0 {
                                    c += nx_i;
                                } else if c >= nx_i {
                                    c -= nx_i;
                                }
                                c as usize
                            };
                            let c = row + cx;
                            let start = offsets[c] as usize;
                            let end = offsets[c + 1] as usize;
                            for cand in &sorted[start..end] {
                                if cand.id == id_u32 {
                                    continue;
                                }
                                let mut dx = cand.x - sx;
                                if dx > half_bx {
                                    dx -= bx;
                                } else if dx < -half_bx {
                                    dx += bx;
                                }
                                let mut dy = cand.y - sy;
                                if dy > half_by {
                                    dy -= by;
                                } else if dy < -half_by {
                                    dy += by;
                                }
                                let d2 = dx * dx + dy * dy;
                                if d2 > min_d2 && d2 <= max_d2 {
                                    buf.push((d2, [dx, dy]));
                                }
                            }
                        }
                    }
                };

                // Stage 1: inner radius (bounds >99% of Voronoi cells completely).
                gather(
                    rx1,
                    ry1,
                    full_x1,
                    full_y1,
                    1e-24,
                    r_inner2,
                    &mut scratch.neighbours,
                );
                scratch
                    .neighbours
                    .sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

                scratch.poly.clear();
                scratch.poly.extend_from_slice(&[
                    [-r_outer, -r_outer],
                    [r_outer, -r_outer],
                    [r_outer, r_outer],
                    [-r_outer, r_outer],
                ]);
                let mut r_max_sq = 2.0 * r_outer2;
                if !clip_sorted_plane(
                    &scratch.neighbours,
                    &mut scratch.poly,
                    &mut scratch.next,
                    &mut r_max_sq,
                ) {
                    *out = [0.0, 0.0];
                    return;
                }

                // Stage 2 (rare): if a site beyond r_inner could still clip the
                // polygon (2 * r_max > r_inner), scan the outer shell.
                if !same_radius && 4.0 * r_max_sq * (1.0 + 1e-12) > r_inner2 {
                    gather(
                        rx2,
                        ry2,
                        full_x2,
                        full_y2,
                        r_inner2,
                        r_outer2,
                        &mut scratch.neighbours,
                    );
                    scratch
                        .neighbours
                        .sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
                    if !clip_sorted_plane(
                        &scratch.neighbours,
                        &mut scratch.poly,
                        &mut scratch.next,
                        &mut r_max_sq,
                    ) {
                        *out = [0.0, 0.0];
                        return;
                    }
                }

                *out = polygon_centroid(&scratch.poly).unwrap_or([0.0, 0.0]);
            });

        let mut max_disp_sq = 0.0_f64;
        for (site, &[dx, dy]) in sites[..n_movable].iter_mut().zip(&self.disps) {
            max_disp_sq = max_disp_sq.max(dx * dx + dy * dy);
            site[0] = wrap(site[0] + dx, bx);
            site[1] = wrap(site[1] + dy, by);
        }
        max_disp_sq.sqrt()
    }

    /// Return the per-site displacements computed by the most recent [`step_sites`].
    pub(crate) fn last_displacements(&self) -> &[[f64; 2]] {
        &self.disps
    }
}

/// One planar Lloyd step for a group of lipids sharing a periodic box `box_xy`
/// (nm): move each lipid rigidly in x/y so its lateral centroid lands on its
/// Voronoi-cell centroid, then wrap into the box. Beads keep their z, so every
/// lipid stays in its leaflet. Returns the largest in-plane displacement (nm).
///
/// Deterministic: the per-site cell computations are independent and written in
/// index order, so a parallel run is bit-identical to a serial one.
pub fn cvt_plane_step(
    positions: &mut [[f64; 3]],
    lipids: &[Range<usize>],
    box_xy: [f64; 2],
) -> f64 {
    if lipids.len() < 2 || box_xy[0] <= 0.0 || box_xy[1] <= 0.0 {
        return 0.0;
    }
    let mut sites = extract_plane_sites(positions, lipids, box_xy);
    let mut ws = PlaneCvtWorkspace::new();
    let max_disp = ws.step_sites(&mut sites, box_xy);
    for (range, &[dx, dy]) in lipids.iter().zip(ws.last_displacements()) {
        for k in range.clone() {
            positions[k][0] = wrap(positions[k][0] + dx, box_xy[0]);
            positions[k][1] = wrap(positions[k][1] + dy, box_xy[1]);
        }
    }
    max_disp
}

// ============================ Spherical ============================

#[inline]
pub(crate) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let m = norm3(a);
    if m < 1e-12 {
        a
    } else {
        let inv = 1.0 / m;
        [a[0] * inv, a[1] * inv, a[2] * inv]
    }
}

/// Rotate `v` about a unit `axis` through the origin (Rodrigues), given
/// precomputed `(sin_a, cos_a)`.
#[inline]
pub(crate) fn rotate_with_sin_cos(v: [f64; 3], axis: [f64; 3], sin_a: f64, cos_a: f64) -> [f64; 3] {
    let cross = cross3(axis, v);
    let dot_k = dot3(axis, v) * (1.0 - cos_a);
    [
        v[0] * cos_a + cross[0] * sin_a + axis[0] * dot_k,
        v[1] * cos_a + cross[1] * sin_a + axis[1] * dot_k,
        v[2] * cos_a + cross[2] * sin_a + axis[2] * dot_k,
    ]
}

/// Rotate every bead in `range` about `center` by the minimal geodesic rotation
/// carrying unit vector `from` to unit vector `to`.
pub(crate) fn rotate_lipid_from_to(
    positions: &mut [[f64; 3]],
    range: Range<usize>,
    center: [f64; 3],
    from: [f64; 3],
    to: [f64; 3],
) {
    let axis = cross3(from, to);
    let sin = norm3(axis);
    let cos = dot3(from, to).clamp(-1.0, 1.0);
    if sin < 1e-12 {
        if cos < 0.0 {
            // 180-degree flip: rotate about any axis perpendicular to `from`.
            let perp = any_perpendicular(from);
            for k in range {
                let rel = [
                    positions[k][0] - center[0],
                    positions[k][1] - center[1],
                    positions[k][2] - center[2],
                ];
                let rot = rotate_with_sin_cos(rel, perp, 0.0, -1.0);
                positions[k] = [center[0] + rot[0], center[1] + rot[1], center[2] + rot[2]];
            }
        }
        return;
    }
    let inv_sin = 1.0 / sin;
    let unit_axis = [axis[0] * inv_sin, axis[1] * inv_sin, axis[2] * inv_sin];
    for k in range {
        let rel = [
            positions[k][0] - center[0],
            positions[k][1] - center[1],
            positions[k][2] - center[2],
        ];
        let rot = rotate_with_sin_cos(rel, unit_axis, sin, cos);
        positions[k] = [center[0] + rot[0], center[1] + rot[1], center[2] + rot[2]];
    }
}

/// A unit vector perpendicular to `v` (used to seed a tangent frame).
#[inline]
fn any_perpendicular(v: [f64; 3]) -> [f64; 3] {
    let other = if v[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    normalize3(cross3(v, other))
}

/// Fill `out` with the 12-gon inscribed in the circle of angular radius `alpha`
/// about unit `site`, using precomputed `(sin_alpha, cos_alpha)` and [`CAP_CS`].
fn fill_cap_polygon(site: [f64; 3], sin_a: f64, cos_a: f64, out: &mut Vec<[f64; 3]>) {
    out.clear();
    let u = any_perpendicular(site);
    let v = cross3(site, u);
    for &(ct, st) in &CAP_CS {
        out.push([
            cos_a * site[0] + sin_a * (ct * u[0] + st * v[0]),
            cos_a * site[1] + sin_a * (ct * u[1] + st * v[1]),
            cos_a * site[2] + sin_a * (ct * u[2] + st * v[2]),
        ]);
    }
}

/// Clip a spherical polygon (unit-vector vertices, great-circle edges) into
/// `out`, keeping `{x : x·normal ≥ 0}`. Each crossed minor arc `a → b` meets the
/// plane at `normalize((1 - t) a + t b)` with `t = sa / (sa - sb)`. Returns the
/// minimum `v·site` (cosine of the angular circumradius) over the output vertices.
fn clip_halfspace_sphere(
    poly: &[[f64; 3]],
    site: [f64; 3],
    normal: [f64; 3],
    out: &mut Vec<[f64; 3]>,
) -> f64 {
    out.clear();
    let n = poly.len();
    if n == 0 {
        return 1.0;
    }
    let mut min_dot = 1.0_f64;
    let mut a = poly[n - 1];
    let mut sa = dot3(a, normal);
    let mut a_in = sa >= 0.0;
    for &b in poly {
        let sb = dot3(b, normal);
        let b_in = sb >= 0.0;
        if a_in != b_in {
            let t = sa / (sa - sb);
            let p = normalize3([
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ]);
            min_dot = min_dot.min(dot3(p, site));
            out.push(p);
        }
        if b_in {
            min_dot = min_dot.min(dot3(b, site));
            out.push(b);
        }
        a = b;
        sa = sb;
        a_in = b_in;
    }
    min_dot
}

/// Clip `poly` against dot-sorted (nearest first) spherical neighbours, stopping
/// once `site·n < cos(2θ_max)`. Updates `min_dot` in place; returns `false` if
/// the cell collapsed below 3 vertices.
fn clip_sorted_sphere(
    site: [f64; 3],
    neighbours: &[(f64, [f64; 3])],
    poly: &mut Vec<[f64; 3]>,
    next: &mut Vec<[f64; 3]>,
    min_dot: &mut f64,
) -> bool {
    for &(dot_sn, n) in neighbours {
        if *min_dot > 0.0 && dot_sn < 2.0 * *min_dot * *min_dot - 1.0 - 1e-12 {
            break;
        }
        // Half-space closer to `site`: x·(site − n) ≥ 0.
        let normal = [site[0] - n[0], site[1] - n[1], site[2] - n[2]];
        *min_dot = clip_halfspace_sphere(poly, site, normal, next);
        std::mem::swap(poly, next);
        if poly.len() < 3 {
            return false;
        }
    }
    true
}

/// Solid angle (area on the unit sphere) of the spherical triangle `(a, b, c)` by
/// the Van Oosterom–Strackee formula.
#[inline]
fn spherical_triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let numerator = dot3(a, cross3(b, c)).abs();
    let denominator = 1.0 + dot3(a, b) + dot3(b, c) + dot3(c, a);
    2.0 * numerator.atan2(denominator)
}

/// The area-weighted centroid (a unit direction) of a spherical polygon, fanned
/// from its first vertex into triangles. `None` for a degenerate cell.
fn spherical_polygon_centroid(poly: &[[f64; 3]]) -> Option<[f64; 3]> {
    if poly.len() < 3 {
        return None;
    }
    let mut acc = [0.0; 3];
    let mut total = 0.0;
    let a = poly[0];
    for i in 1..poly.len() - 1 {
        let (b, c) = (poly[i], poly[i + 1]);
        let area = spherical_triangle_area(a, b, c);
        let centroid = normalize3([a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]]);
        acc[0] += area * centroid[0];
        acc[1] += area * centroid[1];
        acc[2] += area * centroid[2];
        total += area;
    }
    if total < 1e-12 {
        return None;
    }
    Some(normalize3(acc))
}

#[derive(Clone, Copy, Default)]
struct SphereSite {
    dir: [f64; 3],
    id: u32,
}

struct SphereScratch {
    neighbours: Vec<(f64, [f64; 3])>,
    poly: Vec<[f64; 3]>,
    next: Vec<[f64; 3]>,
}

impl SphereScratch {
    fn new() -> Self {
        Self {
            neighbours: Vec::with_capacity(32),
            poly: Vec::with_capacity(20),
            next: Vec::with_capacity(20),
        }
    }
}

/// Extract each lipid's unit centroid direction from `center`.
pub(crate) fn extract_sphere_dirs(
    positions: &[[f64; 3]],
    lipids: &[Range<usize>],
    center: [f64; 3],
) -> Vec<[f64; 3]> {
    lipids
        .iter()
        .map(|range| {
            let m = range.len() as f64;
            let mut c = [0.0; 3];
            for k in range.clone() {
                c[0] += positions[k][0];
                c[1] += positions[k][1];
                c[2] += positions[k][2];
            }
            normalize3([
                c[0] / m - center[0],
                c[1] / m - center[1],
                c[2] / m - center[2],
            ])
        })
        .collect()
}

/// Reusable scratch buffers for multi-round spherical CVT (flat 3D CSR grid over
/// the unit cube `[-1, 1]³` with two-stage adaptive radius search).
pub(crate) struct SphereCvtWorkspace {
    cell_coords: Vec<[u16; 3]>,
    counts: Vec<u32>,
    offsets: Vec<u32>,
    sorted: Vec<SphereSite>,
    targets: Vec<[f64; 3]>,
}

impl SphereCvtWorkspace {
    pub(crate) fn new() -> Self {
        Self {
            cell_coords: Vec::new(),
            counts: Vec::new(),
            offsets: Vec::new(),
            sorted: Vec::new(),
            targets: Vec::new(),
        }
    }

    /// Advance unit shell directions `dirs` by one exact spherical Lloyd step in
    /// place. Returns the largest surface displacement `radius * max_angle` (nm).
    pub(crate) fn step_dirs(&mut self, dirs: &mut [[f64; 3]], radius: f64) -> f64 {
        self.step_dirs_pinned(dirs, dirs.len(), radius)
    }

    /// Like [`step_dirs`](Self::step_dirs), but only `dirs[..n_movable]` move;
    /// `dirs[n_movable..]` stay fixed as obstacle boundary directions that still
    /// clip movable spherical Voronoi cells.
    pub(crate) fn step_dirs_pinned(
        &mut self,
        dirs: &mut [[f64; 3]],
        n_movable: usize,
        radius: f64,
    ) -> f64 {
        let n = dirs.len();
        let n_movable = n_movable.min(n);
        if n_movable == 0 || n < 2 || radius <= 0.0 {
            return 0.0;
        }

        let beta = (4.0 * std::f64::consts::PI / n as f64).sqrt();
        let alpha_outer = (OUTER_SPACINGS * beta).min(MAX_CAP_ALPHA);
        let alpha_inner = (INNER_SPACINGS * beta).min(alpha_outer);
        let (sin_outer, cos_outer) = alpha_outer.sin_cos();
        let cos_inner = alpha_inner.cos();
        let chord_outer = 2.0 * (alpha_outer * 0.5).sin();
        let chord_inner = 2.0 * (alpha_inner * 0.5).sin();

        // Flat 3D grid over [-1, 1]^3 sized to chord_inner (capped at 64^3 cells
        // = 1 MB offsets) so a 3x3x3 scan covers the inner radius completely.
        let dim = ((2.0 / chord_inner.max(1e-6)).floor() as usize).clamp(1, 64);
        let cell_size = 2.0 / dim as f64;
        let inv_cell = dim as f64 * 0.5;
        let max_c = dim - 1;
        let ncells = dim * dim * dim;
        let r1 = ((chord_inner / cell_size).ceil() as usize).max(1);
        let r2 = ((chord_outer / cell_size).ceil() as usize).max(r1);
        let same_radius = alpha_inner >= alpha_outer;

        self.cell_coords.clear();
        self.cell_coords.reserve(n);
        self.counts.resize(ncells, 0);
        self.counts.fill(0);

        for &d in dirs.iter() {
            let cx = (((d[0] + 1.0).max(0.0) * inv_cell) as usize).min(max_c);
            let cy = (((d[1] + 1.0).max(0.0) * inv_cell) as usize).min(max_c);
            let cz = (((d[2] + 1.0).max(0.0) * inv_cell) as usize).min(max_c);
            self.cell_coords.push([cx as u16, cy as u16, cz as u16]);
            self.counts[(cz * dim + cy) * dim + cx] += 1;
        }

        self.offsets.resize(ncells + 1, 0);
        let mut acc = 0u32;
        for c in 0..ncells {
            self.offsets[c] = acc;
            acc += self.counts[c];
            self.counts[c] = self.offsets[c];
        }
        self.offsets[ncells] = acc;

        self.sorted.resize(n, SphereSite::default());
        for i in 0..n {
            let [cx, cy, cz] = self.cell_coords[i];
            let c = (cz as usize * dim + cy as usize) * dim + cx as usize;
            let slot = self.counts[c] as usize;
            self.counts[c] += 1;
            self.sorted[slot] = SphereSite {
                dir: dirs[i],
                id: i as u32,
            };
        }

        let dirs_ref = &*dirs;
        let cell_coords = &self.cell_coords;
        let offsets = &self.offsets;
        let sorted = &self.sorted;

        self.targets.resize(n_movable, [0.0; 3]);
        self.targets
            .par_iter_mut()
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each_init(SphereScratch::new, |scratch, (i, out)| {
                let s = dirs_ref[i];
                let [cx, cy, cz] = cell_coords[i];
                let (cx, cy, cz) = (cx as usize, cy as usize, cz as usize);
                let id_u32 = i as u32;

                let gather = |r_cells: usize,
                              min_dot_cut: f64,
                              max_dot_cut: f64,
                              buf: &mut Vec<(f64, [f64; 3])>| {
                    buf.clear();
                    let x0 = cx.saturating_sub(r_cells);
                    let x1 = (cx + r_cells).min(max_c);
                    let y0 = cy.saturating_sub(r_cells);
                    let y1 = (cy + r_cells).min(max_c);
                    let z0 = cz.saturating_sub(r_cells);
                    let z1 = (cz + r_cells).min(max_c);
                    for z in z0..=z1 {
                        let z_base = z * dim;
                        for y in y0..=y1 {
                            let row = (z_base + y) * dim;
                            // x is the contiguous axis: scan x0..=x1 as one slice.
                            let start = offsets[row + x0] as usize;
                            let end = offsets[row + x1 + 1] as usize;
                            for cand in &sorted[start..end] {
                                if cand.id != id_u32 {
                                    let d = dot3(s, cand.dir);
                                    if d >= min_dot_cut && d < max_dot_cut {
                                        buf.push((d, cand.dir));
                                    }
                                }
                            }
                        }
                    }
                };

                // Stage 1: inner angular radius (bounds >99% of Voronoi cells).
                gather(r1, cos_inner, 2.0, &mut scratch.neighbours);
                scratch
                    .neighbours
                    .sort_unstable_by(|a, b| b.0.total_cmp(&a.0));

                fill_cap_polygon(s, sin_outer, cos_outer, &mut scratch.poly);
                let mut min_dot = cos_outer;
                if !clip_sorted_sphere(
                    s,
                    &scratch.neighbours,
                    &mut scratch.poly,
                    &mut scratch.next,
                    &mut min_dot,
                ) {
                    *out = s;
                    return;
                }

                // Stage 2 (rare): if a site outside alpha_inner could still clip
                // the cell (2 * theta_max > alpha_inner), scan the outer shell.
                if !same_radius
                    && (min_dot <= 0.0 || 2.0 * min_dot * min_dot - 1.0 - 1e-12 < cos_inner)
                {
                    gather(r2, cos_outer, cos_inner, &mut scratch.neighbours);
                    scratch
                        .neighbours
                        .sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
                    if !clip_sorted_sphere(
                        s,
                        &scratch.neighbours,
                        &mut scratch.poly,
                        &mut scratch.next,
                        &mut min_dot,
                    ) {
                        *out = s;
                        return;
                    }
                }

                *out = spherical_polygon_centroid(&scratch.poly).unwrap_or(s);
            });

        // Compute the largest angular move with a single acos() per round.
        let mut min_cos = 1.0_f64;
        for (from, &to) in dirs[..n_movable].iter_mut().zip(&self.targets) {
            min_cos = min_cos.min(dot3(*from, to));
            *from = to;
        }
        radius * min_cos.clamp(-1.0, 1.0).acos()
    }
}

/// One spherical Lloyd step for a group of lipids on a shell of `center`/`radius`:
/// move each lipid rigidly (a rotation about `center`) so its centroid direction
/// lands on its spherical Voronoi-cell centroid, keeping every bead on the shell.
/// Returns the largest surface displacement (nm). Deterministic and parallel.
pub fn cvt_sphere_step(
    positions: &mut [[f64; 3]],
    lipids: &[Range<usize>],
    center: [f64; 3],
    radius: f64,
) -> f64 {
    if lipids.len() < 2 || radius <= 0.0 {
        return 0.0;
    }
    let init_dirs = extract_sphere_dirs(positions, lipids, center);
    let mut dirs = init_dirs.clone();
    let mut ws = SphereCvtWorkspace::new();
    let max_disp = ws.step_dirs(&mut dirs, radius);
    for (range, (&from, &to)) in lipids.iter().zip(init_dirs.iter().zip(&dirs)) {
        rotate_lipid_from_to(positions, range.clone(), center, from, to);
    }
    max_disp
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell_centroid_plane(site: [f64; 2], neighbours: &[[f64; 2]], half: f64) -> [f64; 2] {
        let mut rel: Vec<(f64, [f64; 2])> = neighbours
            .iter()
            .map(|&n| {
                let d = [n[0] - site[0], n[1] - site[1]];
                (d[0] * d[0] + d[1] * d[1], d)
            })
            .collect();
        rel.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let mut poly = vec![[-half, -half], [half, -half], [half, half], [-half, half]];
        let mut next = Vec::new();
        let mut r_max_sq = 2.0 * half * half;
        if !clip_sorted_plane(&rel, &mut poly, &mut next, &mut r_max_sq) {
            return site;
        }
        let disp = polygon_centroid(&poly).unwrap_or([0.0, 0.0]);
        [site[0] + disp[0], site[1] + disp[1]]
    }

    fn cell_centroid_sphere(site: [f64; 3], neighbours: &[[f64; 3]], alpha: f64) -> [f64; 3] {
        let mut sorted: Vec<(f64, [f64; 3])> =
            neighbours.iter().map(|&n| (dot3(site, n), n)).collect();
        sorted.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
        let (sa, ca) = alpha.sin_cos();
        let mut poly = Vec::new();
        let mut next = Vec::new();
        fill_cap_polygon(site, sa, ca, &mut poly);
        let mut min_dot = ca;
        if !clip_sorted_sphere(site, &sorted, &mut poly, &mut next, &mut min_dot) {
            return site;
        }
        spherical_polygon_centroid(&poly).unwrap_or(site)
    }

    #[test]
    fn clip_to_a_unit_square() {
        // Four neighbours two units away on the axes: the cell is the unit square
        // centred on the origin, so its centroid is the origin.
        let site = [0.0, 0.0];
        let neighbours = [[2.0, 0.0], [-2.0, 0.0], [0.0, 2.0], [0.0, -2.0]];
        let c = cell_centroid_plane(site, &neighbours, 10.0);
        assert!(c[0].abs() < 1e-9 && c[1].abs() < 1e-9, "centroid {c:?}");
    }

    #[test]
    fn centroid_pulls_toward_open_space() {
        // Neighbours only on the +x side: the cell opens toward −x, so the centroid
        // moves to negative x.
        let site = [0.0, 0.0];
        let neighbours = [[2.0, 0.0], [1.6, 1.6], [1.6, -1.6]];
        let c = cell_centroid_plane(site, &neighbours, 6.0);
        assert!(c[0] < -0.05, "expected a pull toward −x, got {c:?}");
    }

    #[test]
    fn polygon_centroid_of_a_unit_square() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let c = polygon_centroid(&square).unwrap();
        assert!((c[0] - 0.5).abs() < 1e-12 && (c[1] - 0.5).abs() < 1e-12);
    }

    /// A perfect periodic lattice is already centroidal: a Lloyd step barely moves.
    #[test]
    fn uniform_lattice_is_a_fixed_point() {
        let m = 12;
        let step = 1.0;
        let box_xy = [m as f64 * step, m as f64 * step];
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        for i in 0..m {
            for j in 0..m {
                let k = positions.len();
                positions.push([i as f64 * step, j as f64 * step, 0.0]);
                lipids.push(k..k + 1);
            }
        }
        let moved = cvt_plane_step(&mut positions, &lipids, box_xy);
        assert!(moved < 1e-6 * step, "lattice drifted by {moved}");
    }

    /// Lloyd steps relax a clustered set toward uniform: the mean nearest-neighbour
    /// distance grows toward the ideal lattice spacing as the cluster spreads.
    #[test]
    fn clustered_points_spread_out() {
        let box_xy = [16.0, 16.0];
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        // 64 lipids bunched into a 4×4 corner block (spacing 0.5), far from uniform.
        for i in 0..8 {
            for j in 0..8 {
                let k = positions.len();
                positions.push([1.0 + i as f64 * 0.5, 1.0 + j as f64 * 0.5, 0.0]);
                lipids.push(k..k + 1);
            }
        }
        let nn_mean = |pos: &[[f64; 3]]| {
            let mut total = 0.0;
            for a in 0..pos.len() {
                let mut best = f64::MAX;
                for b in 0..pos.len() {
                    if a == b {
                        continue;
                    }
                    let mut dx = pos[a][0] - pos[b][0];
                    let mut dy = pos[a][1] - pos[b][1];
                    dx -= box_xy[0] * (dx / box_xy[0]).round();
                    dy -= box_xy[1] * (dy / box_xy[1]).round();
                    best = best.min(dx * dx + dy * dy);
                }
                total += best.sqrt();
            }
            total / pos.len() as f64
        };
        let before = nn_mean(&positions);
        for _ in 0..40 {
            cvt_plane_step(&mut positions, &lipids, box_xy);
        }
        let after = nn_mean(&positions);
        // 64 lipids in a 16×16 box ⇒ ideal spacing 2.0; the cluster must open up.
        assert!(
            after > before + 0.5,
            "cluster did not spread: nn {before:.2} -> {after:.2}"
        );
    }

    /// Same seed of positions ⇒ identical relaxation (parallel == serial order).
    #[test]
    fn cvt_step_is_deterministic() {
        let box_xy = [10.0, 10.0];
        let make = || {
            let mut positions = Vec::new();
            let mut lipids = Vec::new();
            // A jittered grid (deterministic, no RNG): reproducible across runs.
            for i in 0..10 {
                for j in 0..10 {
                    let k = positions.len();
                    let jx = ((i * 7 + j * 3) % 5) as f64 * 0.05;
                    let jy = ((i * 3 + j * 5) % 5) as f64 * 0.05;
                    positions.push([i as f64 + jx, j as f64 + jy, 0.0]);
                    lipids.push(k..k + 1);
                }
            }
            (positions, lipids)
        };
        let (mut a, la) = make();
        let (mut b, lb) = make();
        for _ in 0..5 {
            cvt_plane_step(&mut a, &la, box_xy);
            cvt_plane_step(&mut b, &lb, box_xy);
        }
        assert_eq!(a, b);
    }

    // ---- spherical ----

    #[test]
    fn octahedron_cell_centroid_is_the_site() {
        // +x with the four equatorial neighbours ±y, ±z: by symmetry its spherical
        // Voronoi cell is centred on +x, so the centroid returns to the site.
        let site = [1.0, 0.0, 0.0];
        let neighbours = [
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let c = cell_centroid_sphere(site, &neighbours, 1.3);
        assert!((c[0] - 1.0).abs() < 1e-6, "centroid {c:?} not on +x");
        assert!(
            c[1].abs() < 1e-6 && c[2].abs() < 1e-6,
            "centroid {c:?} off-axis"
        );
    }

    /// Single-bead lipids on a shell, clustered near +z, with the radius `r`.
    #[allow(clippy::type_complexity)]
    fn shell_cluster() -> (Vec<[f64; 3]>, Vec<Range<usize>>, [f64; 3], f64) {
        let center = [0.0, 0.0, 0.0];
        let r = 10.0;
        let mut positions = Vec::new();
        let mut lipids = Vec::new();
        for i in 0..6 {
            for j in 0..7 {
                let a = (i as f64 - 2.5) * 0.08;
                let b = (j as f64 - 3.0) * 0.08;
                let d = normalize3([a, b, 1.0]); // a ~0.3 rad cap about +z
                let k = positions.len();
                positions.push([d[0] * r, d[1] * r, d[2] * r]);
                lipids.push(k..k + 1);
            }
        }
        (positions, lipids, center, r)
    }

    #[test]
    fn spherical_cluster_spreads_over_the_sphere() {
        let (mut positions, lipids, center, r) = shell_cluster();
        // Concentration = |mean direction|: 1 when bunched, →0 as it covers the sphere.
        let concentration = |pos: &[[f64; 3]]| {
            let mut s = [0.0; 3];
            for p in pos {
                let d = normalize3([p[0] - center[0], p[1] - center[1], p[2] - center[2]]);
                for a in 0..3 {
                    s[a] += d[a];
                }
            }
            norm3(s) / pos.len() as f64
        };
        let before = concentration(&positions);
        for _ in 0..40 {
            cvt_sphere_step(&mut positions, &lipids, center, r);
        }
        let after = concentration(&positions);
        assert!(
            after < before - 0.2,
            "cluster did not spread: concentration {before:.2} -> {after:.2}"
        );
    }

    #[test]
    fn spherical_step_keeps_beads_on_the_shell() {
        let (mut positions, lipids, center, r) = shell_cluster();
        cvt_sphere_step(&mut positions, &lipids, center, r);
        for p in &positions {
            let radius = norm3([p[0] - center[0], p[1] - center[1], p[2] - center[2]]);
            assert!(
                (radius - r).abs() < 1e-9,
                "bead left the shell: r = {radius}"
            );
        }
    }
}
