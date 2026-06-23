//! Centroidal Voronoi tessellation (Lloyd's relaxation) by **per-site clipping**.
//!
//! A site's Voronoi cell is, by definition, the intersection of the half-spaces
//! `{x : x is closer to the site than to neighbour n}` over all neighbours `n`.
//! Rather than build a global Delaunay/convex-hull triangulation, each cell is
//! computed independently by clipping a starting polygon against the perpendicular
//! bisectors of nearby sites (Sutherland–Hodgman). This yields the exact Voronoi
//! cell, is robust (no global topology to get wrong), parallelises per site, and
//! needs no external geometry dependency.
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

// ============================ Planar ============================

/// Keep the part of convex polygon `poly` on the *site* side of the perpendicular
/// bisector between the site and a neighbour: the half-plane `{p : (p − m)·dir ≤ 0}`
/// where `m` is the bisector midpoint and `dir` points from the site to the
/// neighbour (so the site itself, at `(site − m)·dir < 0`, is kept).
fn clip_halfplane(poly: &[[f64; 2]], m: [f64; 2], dir: [f64; 2]) -> Vec<[f64; 2]> {
    if poly.is_empty() {
        return Vec::new();
    }
    let signed = |p: [f64; 2]| (p[0] - m[0]) * dir[0] + (p[1] - m[1]) * dir[1];
    let mut out = Vec::with_capacity(poly.len() + 1);
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        let (sa, sb) = (signed(a), signed(b));
        let a_in = sa <= 0.0;
        let b_in = sb <= 0.0;
        if a_in {
            out.push(a);
        }
        // Crossing the boundary: add the intersection point.
        if a_in != b_in {
            let t = sa / (sa - sb);
            out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
    }
    out
}

/// Area centroid of a simple polygon (the shoelace formula). `None` when the
/// polygon is degenerate (near-zero area).
fn polygon_centroid(poly: &[[f64; 2]]) -> Option<[f64; 2]> {
    if poly.len() < 3 {
        return None;
    }
    let (mut a2, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for i in 0..poly.len() {
        let p = poly[i];
        let q = poly[(i + 1) % poly.len()];
        let cross = p[0] * q[1] - q[0] * p[1];
        a2 += cross;
        cx += (p[0] + q[0]) * cross;
        cy += (p[1] + q[1]) * cross;
    }
    if a2.abs() < 1e-12 {
        return None;
    }
    Some([cx / (3.0 * a2), cy / (3.0 * a2)])
}

/// The Voronoi-cell centroid of `site` given `neighbours` (their positions, already
/// expressed relative to `site` under the minimum-image convention). The cell is
/// the `half`-sided starting square around `site` clipped by every neighbour's
/// bisector. Falls back to `site` if the cell collapses.
fn cell_centroid_plane(site: [f64; 2], neighbours: &[[f64; 2]], half: f64) -> [f64; 2] {
    let mut poly = vec![
        [site[0] - half, site[1] - half],
        [site[0] + half, site[1] - half],
        [site[0] + half, site[1] + half],
        [site[0] - half, site[1] + half],
    ];
    for &n in neighbours {
        let m = [(site[0] + n[0]) * 0.5, (site[1] + n[1]) * 0.5];
        let dir = [n[0] - site[0], n[1] - site[1]];
        poly = clip_halfplane(&poly, m, dir);
        if poly.len() < 3 {
            return site;
        }
    }
    polygon_centroid(&poly).unwrap_or(site)
}

/// A uniform cell grid over a periodic box for radius-bounded neighbour queries.
struct CellGridPlane {
    box_xy: [f64; 2],
    nx: usize,
    ny: usize,
    cell_x: f64,
    cell_y: f64,
    buckets: Vec<Vec<usize>>,
}

impl CellGridPlane {
    /// Bin `points` into cells about `target` wide (at least one cell per axis).
    fn new(points: &[[f64; 2]], box_xy: [f64; 2], target: f64) -> Self {
        let nx = ((box_xy[0] / target).floor() as usize).max(1);
        let ny = ((box_xy[1] / target).floor() as usize).max(1);
        let cell_x = box_xy[0] / nx as f64;
        let cell_y = box_xy[1] / ny as f64;
        let mut buckets = vec![Vec::new(); nx * ny];
        for (i, p) in points.iter().enumerate() {
            let gx = ((p[0].rem_euclid(box_xy[0])) / cell_x).floor() as usize % nx;
            let gy = ((p[1].rem_euclid(box_xy[1])) / cell_y).floor() as usize % ny;
            buckets[gy * nx + gx].push(i);
        }
        Self {
            box_xy,
            nx,
            ny,
            cell_x,
            cell_y,
            buckets,
        }
    }

    /// Sites within `radius` of `points[i]`, as positions relative to it under the
    /// minimum-image convention (so the caller sees one image per neighbour;
    /// `radius` must be ≤ half the box so that image is the only one in range).
    fn neighbours(&self, points: &[[f64; 2]], i: usize, radius: f64) -> Vec<[f64; 2]> {
        let s = points[i];
        let gx = ((s[0].rem_euclid(self.box_xy[0])) / self.cell_x).floor() as isize;
        let gy = ((s[1].rem_euclid(self.box_xy[1])) / self.cell_y).floor() as isize;
        let rx = (radius / self.cell_x).ceil() as isize;
        let ry = (radius / self.cell_y).ceil() as isize;
        let r2 = radius * radius;
        let mut out = Vec::new();
        for dy in -ry..=ry {
            for dx in -rx..=rx {
                let cx = (gx + dx).rem_euclid(self.nx as isize) as usize;
                let cy = (gy + dy).rem_euclid(self.ny as isize) as usize;
                for &j in &self.buckets[cy * self.nx + cx] {
                    if j == i {
                        continue;
                    }
                    let mut d = [points[j][0] - s[0], points[j][1] - s[1]];
                    d[0] -= self.box_xy[0] * (d[0] / self.box_xy[0]).round();
                    d[1] -= self.box_xy[1] * (d[1] / self.box_xy[1]).round();
                    if d[0] * d[0] + d[1] * d[1] <= r2 {
                        out.push([s[0] + d[0], s[1] + d[1]]);
                    }
                }
            }
        }
        out
    }
}

/// The lateral (x, y) centroid of lipid `range`'s beads.
fn lipid_centroid_xy(positions: &[[f64; 3]], range: &Range<usize>) -> [f64; 2] {
    let n = range.len() as f64;
    let (mut x, mut y) = (0.0, 0.0);
    for i in range.clone() {
        x += positions[i][0];
        y += positions[i][1];
    }
    [x / n, y / n]
}

/// One planar Lloyd step for a group of lipids sharing a periodic box `box_xy`
/// (nm): move each lipid rigidly in x/y so its lateral centroid lands on its
/// Voronoi-cell centroid, then wrap into the box. Beads keep their z, so every
/// lipid stays in its leaflet. Returns the largest in-plane displacement (nm).
///
/// Deterministic: the per-site cell computations are independent and collected in
/// order, so a parallel run is bit-identical to a serial one.
pub fn cvt_plane_step(
    positions: &mut [[f64; 3]],
    lipids: &[Range<usize>],
    box_xy: [f64; 2],
) -> f64 {
    let n = lipids.len();
    if n < 2 || box_xy[0] <= 0.0 || box_xy[1] <= 0.0 {
        return 0.0;
    }
    let sites: Vec<[f64; 2]> = lipids
        .iter()
        .map(|r| lipid_centroid_xy(positions, r))
        .collect();
    let spacing = (box_xy[0] * box_xy[1] / n as f64).sqrt();
    // Gather neighbours out to a few spacings (enough to bound a near-uniform
    // cell), capped at half the box so the minimum image is unambiguous.
    let radius = (NEIGHBOUR_SPACINGS * spacing).min(0.49 * box_xy[0].min(box_xy[1]));
    let grid = CellGridPlane::new(&sites, box_xy, radius.max(spacing));

    let targets: Vec<[f64; 2]> = (0..n)
        .into_par_iter()
        .map(|i| {
            let neighbours = grid.neighbours(&sites, i, radius);
            cell_centroid_plane(sites[i], &neighbours, radius)
        })
        .collect();

    let mut max_disp = 0.0_f64;
    for (range, (site, target)) in lipids.iter().zip(sites.iter().zip(&targets)) {
        let dx = target[0] - site[0];
        let dy = target[1] - site[1];
        max_disp = max_disp.max((dx * dx + dy * dy).sqrt());
        for k in range.clone() {
            positions[k][0] = (positions[k][0] + dx).rem_euclid(box_xy[0]);
            positions[k][1] = (positions[k][1] + dy).rem_euclid(box_xy[1]);
        }
    }
    max_disp
}

/// Neighbour-search radius, in mean lipid spacings: large enough to bound a cell in
/// a roughly uniform set while keeping the per-site work small.
const NEIGHBOUR_SPACINGS: f64 = 4.0;

// ============================ Spherical ============================

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let m = norm3(a);
    if m < 1e-12 {
        a
    } else {
        [a[0] / m, a[1] / m, a[2] / m]
    }
}

/// Rotate `v` by `angle` (rad) about a unit `axis` through the origin (Rodrigues).
pub(crate) fn rotate_about_axis(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = (angle.sin(), angle.cos());
    let cross = cross3(axis, v);
    let dot = dot3(axis, v);
    [
        v[0] * c + cross[0] * s + axis[0] * dot * (1.0 - c),
        v[1] * c + cross[1] * s + axis[1] * dot * (1.0 - c),
        v[2] * c + cross[2] * s + axis[2] * dot * (1.0 - c),
    ]
}

/// A unit vector perpendicular to `v` (used to seed a tangent frame).
fn any_perpendicular(v: [f64; 3]) -> [f64; 3] {
    let other = if v[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    normalize3(cross3(v, other))
}

/// A geodesic `k`-gon inscribed in the circle of angular radius `alpha` about unit
/// `site`: the starting spherical polygon that the bisectors clip down to the cell.
fn cap_polygon(site: [f64; 3], alpha: f64, k: usize) -> Vec<[f64; 3]> {
    let u = any_perpendicular(site);
    let v = cross3(site, u);
    let (sa, ca) = (alpha.sin(), alpha.cos());
    (0..k)
        .map(|i| {
            let theta = std::f64::consts::TAU * i as f64 / k as f64;
            let (st, ct) = (theta.sin(), theta.cos());
            normalize3([
                ca * site[0] + sa * (ct * u[0] + st * v[0]),
                ca * site[1] + sa * (ct * u[1] + st * v[1]),
                ca * site[2] + sa * (ct * u[2] + st * v[2]),
            ])
        })
        .collect()
}

/// Clip a spherical polygon (unit-vector vertices, great-circle edges) by the
/// half-space `{x : x·normal ≥ 0}` (a plane through the origin). Sutherland–Hodgman,
/// with each crossed edge cut at the great-circle/plane intersection.
fn clip_halfspace_sphere(poly: &[[f64; 3]], normal: [f64; 3]) -> Vec<[f64; 3]> {
    if poly.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(poly.len() + 1);
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        let (sa, sb) = (dot3(a, normal), dot3(b, normal));
        let a_in = sa >= 0.0;
        let b_in = sb >= 0.0;
        if a_in {
            out.push(a);
        }
        if a_in != b_in {
            // Intersection of the edge's great circle (plane normal a×b) with the
            // clip plane (normal `normal`): the line common to both planes.
            let g = cross3(a, b);
            let mut p = normalize3(cross3(g, normal));
            // Pick the root on the minor arc a→b (the one toward the arc midpoint).
            if dot3(p, [a[0] + b[0], a[1] + b[1], a[2] + b[2]]) < 0.0 {
                p = [-p[0], -p[1], -p[2]];
            }
            out.push(p);
        }
    }
    out
}

/// Solid angle (area on the unit sphere) of the spherical triangle `(a, b, c)` by
/// the Van Oosterom–Strackee formula.
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
    for i in 1..poly.len() - 1 {
        let (a, b, c) = (poly[0], poly[i], poly[i + 1]);
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

/// The Voronoi-cell centroid direction of unit `site` given unit `neighbours`,
/// computed by clipping the cap of angular radius `alpha` about `site` by every
/// neighbour's great-circle bisector. Falls back to `site` if the cell collapses.
fn cell_centroid_sphere(site: [f64; 3], neighbours: &[[f64; 3]], alpha: f64) -> [f64; 3] {
    let mut poly = cap_polygon(site, alpha, CAP_VERTICES);
    for &n in neighbours {
        // Half-space closer to `site`: x·(site − n) ≥ 0.
        let normal = [site[0] - n[0], site[1] - n[1], site[2] - n[2]];
        poly = clip_halfspace_sphere(&poly, normal);
        if poly.len() < 3 {
            return site;
        }
    }
    spherical_polygon_centroid(&poly).unwrap_or(site)
}

/// Vertices of the starting cap polygon: enough that its inradius comfortably
/// exceeds a cell, so a cell never pokes through the cap before clipping.
const CAP_VERTICES: usize = 12;

/// A sparse 3D cell grid over unit direction vectors for angular neighbour queries.
struct CellGridSphere {
    cell: f64,
    buckets: std::collections::HashMap<[i64; 3], Vec<usize>>,
}

impl CellGridSphere {
    /// Bin unit `dirs` into cells of side `cell` (a chord length).
    fn new(dirs: &[[f64; 3]], cell: f64) -> Self {
        let cell = cell.max(1e-6);
        let mut buckets: std::collections::HashMap<[i64; 3], Vec<usize>> =
            std::collections::HashMap::new();
        for (i, d) in dirs.iter().enumerate() {
            buckets.entry(Self::key(*d, cell)).or_default().push(i);
        }
        Self { cell, buckets }
    }

    fn key(d: [f64; 3], cell: f64) -> [i64; 3] {
        [
            (d[0] / cell).floor() as i64,
            (d[1] / cell).floor() as i64,
            (d[2] / cell).floor() as i64,
        ]
    }

    /// Directions within angular `alpha` of `dirs[i]` (i.e. dot ≥ cos α), excluding
    /// `i`. Scans the 3×3×3 cell block, which spans one chord ≥ the search radius.
    fn neighbours(&self, dirs: &[[f64; 3]], i: usize, alpha: f64) -> Vec<[f64; 3]> {
        let s = dirs[i];
        let key = Self::key(s, self.cell);
        let cos_alpha = alpha.cos();
        let mut out = Vec::new();
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let k = [key[0] + dx, key[1] + dy, key[2] + dz];
                    if let Some(bucket) = self.buckets.get(&k) {
                        for &j in bucket {
                            if j != i && dot3(s, dirs[j]) >= cos_alpha {
                                out.push(dirs[j]);
                            }
                        }
                    }
                }
            }
        }
        out
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
    let n = lipids.len();
    if n < 2 || radius <= 0.0 {
        return 0.0;
    }
    // Each lipid's centroid direction from the sphere centre.
    let dirs: Vec<[f64; 3]> = lipids
        .iter()
        .map(|r| {
            let m = r.len() as f64;
            let mut c = [0.0; 3];
            for k in r.clone() {
                for a in 0..3 {
                    c[a] += positions[k][a];
                }
            }
            normalize3([
                c[0] / m - center[0],
                c[1] / m - center[1],
                c[2] / m - center[2],
            ])
        })
        .collect();

    // Angular spacing on the shell, and the neighbour/cap angular radius.
    let beta = (4.0 * std::f64::consts::PI / n as f64).sqrt();
    let alpha = (NEIGHBOUR_SPACINGS * beta).min(MAX_CAP_ALPHA);
    let chord = 2.0 * (alpha * 0.5).sin();
    let grid = CellGridSphere::new(&dirs, chord);

    let targets: Vec<[f64; 3]> = (0..n)
        .into_par_iter()
        .map(|i| {
            let neighbours = grid.neighbours(&dirs, i, alpha);
            cell_centroid_sphere(dirs[i], &neighbours, alpha)
        })
        .collect();

    let mut max_disp = 0.0_f64;
    for (range, (from, to)) in lipids.iter().zip(dirs.iter().zip(&targets)) {
        let axis = cross3(*from, *to);
        let sin = norm3(axis);
        let cos = dot3(*from, *to).clamp(-1.0, 1.0);
        let angle = sin.atan2(cos);
        if sin < 1e-12 || angle < 1e-9 {
            continue; // already centroidal (or degenerate)
        }
        max_disp = max_disp.max(radius * angle);
        let axis = [axis[0] / sin, axis[1] / sin, axis[2] / sin];
        for k in range.clone() {
            let rel = [
                positions[k][0] - center[0],
                positions[k][1] - center[1],
                positions[k][2] - center[2],
            ];
            let rot = rotate_about_axis(rel, axis, angle);
            positions[k] = [center[0] + rot[0], center[1] + rot[1], center[2] + rot[2]];
        }
    }
    max_disp
}

/// Largest cap/neighbour angular radius (rad), so a small lipid count (huge `beta`)
/// can't ask for a cap near or past a hemisphere.
const MAX_CAP_ALPHA: f64 = 1.0;

#[cfg(test)]
mod tests {
    use super::*;

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
