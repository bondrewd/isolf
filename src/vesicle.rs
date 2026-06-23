//! Construction of a coarse-grained lipid vesicle (a spherical bilayer).
//!
//! This is the flat-membrane algorithm with the plane replaced by a sphere. Lipid
//! heads are placed at uniform-random directions over an inner and an outer shell,
//! then relaxed to a uniform, clash-free layout by a spherical centroidal-Voronoi
//! (Lloyd) pass per shell ([`crate::pack::equalize_on_sphere`]). Each lipid's beads
//! are stacked **radially** by walking its bond lengths from the head: the outer
//! leaflet points outward (head on the outer surface, tail toward the centre), the
//! inner leaflet points inward, so the two tail regions meet near the mid-radius.
//!
//! Spacing, packing mode, and explicit area overrides are shared with
//! [`crate::membrane`] (the outer leaflet maps to the membrane's "upper",
//! the inner to "lower"), so a vesicle starts at the same density a flat
//! membrane would and is built clash-free.

use std::f64::consts::PI;

use rand::Rng;

use crate::composition::Composition;
use crate::error::BuildError;
use crate::force_field::ForceField;
use crate::membrane::{
    BuildOptions, DEFAULT_PACK_ITERATIONS, DEFAULT_PACK_TOLERANCE, Membrane, Particle, RelaxFrame,
    aggregate_counts, bead_packing_radius, canonicalize, clash_floor, leaflet_spacing,
    molecule_ranges,
};
use crate::pack;

/// The vesicle radius the user fixed. With one radius the other is derived from the
/// bilayer thickness; with both, the gap between them must hold the bilayer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VesicleRadius {
    /// Inner (lumen) radius in nm.
    Inner(f64),
    /// Outer radius in nm.
    Outer(f64),
    /// Both the inner (lumen) and outer radii in nm.
    Both { inner: f64, outer: f64 },
}

/// Build a vesicle from outer and inner leaflet compositions and a radius spec.
///
/// `rng` seeds the random placement; pass a seeded generator for reproducible
/// output.
///
/// # Errors
///
/// Returns a [`BuildError`] if a composition names an unknown lipid or bead, a
/// lipid has no phosphate bead, or the radii are too small to hold the bilayer
/// ([`VesicleTooSmall`](BuildError::VesicleTooSmall) or
/// [`VesicleTooThin`](BuildError::VesicleTooThin)).
pub fn build_vesicle(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: VesicleRadius,
    options: &BuildOptions,
    rng: &mut impl Rng,
) -> Result<Membrane, BuildError> {
    build_vesicle_inner(force_field, name, outer, inner, radius, options, rng, None)
}

/// Like [`build_vesicle`], but also records every relaxation step (the per-shell
/// spherical Lloyd rounds, then the de-clash) as a [`RelaxFrame`], for the `--gif`
/// animation. Each lipid is projected to 2D by a Lambert azimuthal equal-area map,
/// so a uniform spread over the sphere stays uniform in the frame.
///
/// # Errors
///
/// The same as [`build_vesicle`].
pub fn build_vesicle_recorded(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: VesicleRadius,
    options: &BuildOptions,
    rng: &mut impl Rng,
) -> Result<(Membrane, Vec<RelaxFrame>), BuildError> {
    let mut frames = Vec::new();
    let vesicle = build_vesicle_inner(
        force_field,
        name,
        outer,
        inner,
        radius,
        options,
        rng,
        Some(&mut frames),
    )?;
    Ok((vesicle, frames))
}

#[allow(clippy::too_many_arguments)]
fn build_vesicle_inner(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: VesicleRadius,
    options: &BuildOptions,
    rng: &mut impl Rng,
    frames: Option<&mut Vec<RelaxFrame>>,
) -> Result<Membrane, BuildError> {
    let outer_spacing = leaflet_spacing(force_field, outer)?;
    let inner_spacing = leaflet_spacing(force_field, inner)?;

    // Bilayer thickness: each leaflet's longest lipid plus a clash-free gap so
    // the two tail regions do not overlap at the mid-radius.
    let gap = clash_floor(force_field, outer)?.max(clash_floor(force_field, inner)?);
    let thickness =
        max_chain_length(force_field, outer)? + max_chain_length(force_field, inner)? + gap;

    let (inner_radius, outer_radius) = match radius {
        VesicleRadius::Inner(r) => (r, r + thickness),
        VesicleRadius::Outer(r) => {
            if r <= thickness {
                return Err(BuildError::VesicleTooSmall {
                    outer: r,
                    thickness,
                });
            }
            (r - thickness, r)
        }
        VesicleRadius::Both { inner, outer } => {
            if outer - inner < thickness {
                return Err(BuildError::VesicleTooThin {
                    inner,
                    outer,
                    minimum: thickness,
                });
            }
            (inner, outer)
        }
    };

    // Lipids per shell: the shell area divided by the per-lipid area (spacing²).
    let outer_count = shell_count(outer_radius, outer_spacing);
    let inner_count = shell_count(inner_radius, inner_spacing);

    let mut molecule_id = 0;
    let mut atom_id = 0;
    let mut particles = Vec::new();

    let outer_counts = place_shell(
        force_field,
        outer,
        outer_radius,
        outer_count,
        Shell::Outer,
        rng,
        &mut molecule_id,
        &mut atom_id,
        &mut particles,
    )?;
    // The lipids placed so far form the outer shell (the inner follows).
    let outer_lipids = molecule_id;
    let inner_counts = place_shell(
        force_field,
        inner,
        inner_radius,
        inner_count,
        Shell::Inner,
        rng,
        &mut molecule_id,
        &mut atom_id,
        &mut particles,
    )?;

    // Relax the random placement to a uniform, clash-free start: an independent
    // spherical centroidal-Voronoi (Lloyd) pass per shell, then a shared de-clash.
    relax_vesicle(
        force_field,
        &mut particles,
        outer_lipids,
        outer_radius,
        inner_radius,
        frames,
    );

    let lipid_counts = aggregate_counts(outer, &outer_counts, inner, &inner_counts);
    let mut particles = canonicalize(particles, &lipid_counts);

    // Centre the vesicle in a cubic box with vacuum padding on every side; the
    // padding sits well past the force-field cutoff, so periodic images stay apart.
    let box_side = 2.0 * outer_radius + 2.0 * options.padding;
    let centre = 0.5 * box_side;
    for particle in &mut particles {
        for coordinate in &mut particle.position {
            *coordinate += centre;
        }
    }

    Ok(Membrane {
        name: name.to_string(),
        box_size: [box_side; 3],
        lipid_counts,
        particles,
    })
}

/// Which shell a lipid belongs to: the outer one points outward (beads stack
/// toward the centre), the inner one points inward.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Shell {
    Outer,
    Inner,
}

/// Number of lipids that cover a shell of the given radius at `spacing`.
fn shell_count(radius: f64, spacing: f64) -> usize {
    (4.0 * PI * radius * radius / (spacing * spacing)).round() as usize
}

/// The longest lipid chain (sum of bond lengths, nm) in a composition.
fn max_chain_length(
    force_field: &ForceField,
    composition: &Composition,
) -> Result<f64, BuildError> {
    let mut longest = 0.0_f64;
    for component in composition.components() {
        let lipid = force_field
            .lipid(&component.name)
            .ok_or_else(|| BuildError::UnknownLipid(component.name.clone()))?;
        let length: f64 = lipid.bonds.iter().map(|bond| bond.length).sum();
        longest = longest.max(length);
    }
    Ok(longest)
}

/// Place `count` lipids of `composition` over a shell, appending them to
/// `particles` and returning the per-component molecule counts.
#[allow(clippy::too_many_arguments)]
fn place_shell(
    force_field: &ForceField,
    composition: &Composition,
    radius: f64,
    count: usize,
    shell: Shell,
    rng: &mut impl Rng,
    molecule_id: &mut usize,
    atom_id: &mut usize,
    particles: &mut Vec<Particle>,
) -> Result<Vec<usize>, BuildError> {
    let counts = composition.partition(count);
    if count == 0 {
        return Ok(counts);
    }

    // One species per lipid, in partition order; the random directions below mix
    // the species over the sphere, and Lloyd's relaxation makes the layout uniform.
    let mut species: Vec<&str> = Vec::with_capacity(count);
    for (component, &n) in composition.components().iter().zip(&counts) {
        species.extend(std::iter::repeat_n(component.name.as_str(), n));
    }

    for &lipid_name in &species {
        *molecule_id += 1;
        let lipid = force_field
            .lipid(lipid_name)
            .ok_or_else(|| BuildError::UnknownLipid(lipid_name.to_string()))?;
        let direction = random_sphere_direction(rng);

        // Head bead 0 sits on the shell; walk the chain radially toward the
        // bilayer centre (inward for the outer shell, outward for the inner).
        let mut r = radius;
        let last = lipid.beads.len() - 1;
        for (k, bead) in lipid.beads.iter().enumerate() {
            let bead_type = force_field
                .bead_type(bead)
                .ok_or_else(|| BuildError::UnknownBead(bead.clone()))?;
            *atom_id += 1;
            particles.push(Particle {
                residue: lipid.name.clone(),
                bead: bead.clone(),
                molecule_id: *molecule_id,
                atom_id: *atom_id,
                position: [direction[0] * r, direction[1] * r, direction[2] * r],
                charge: bead_type.charge,
                mass: bead_type.mass,
            });
            if k < last {
                r += match shell {
                    Shell::Outer => -lipid.bonds[k].length,
                    Shell::Inner => lipid.bonds[k].length,
                };
            }
        }
    }

    Ok(counts)
}

/// A uniform-random unit direction on the sphere (Archimedes' method: `z` uniform
/// in [−1, 1], azimuth uniform), seeded for reproducibility.
fn random_sphere_direction(rng: &mut impl Rng) -> [f64; 3] {
    let z: f64 = rng.random_range(-1.0..1.0);
    let phi: f64 = rng.random_range(0.0..2.0 * PI);
    let r = (1.0 - z * z).sqrt();
    [r * phi.cos(), r * phi.sin(), z]
}

/// Relax the randomly-placed shells to a uniform, clash-free configuration: an
/// independent spherical centroidal-Voronoi (Lloyd) pass per shell plus a shared
/// de-clash. `outer_lipids` is the number of lipids on the outer shell (placed
/// first), so the two shells' molecule ranges can be separated. Placement is about
/// the origin, so the sphere centre is `[0, 0, 0]`.
fn relax_vesicle(
    force_field: &ForceField,
    particles: &mut [Particle],
    outer_lipids: usize,
    outer_radius: f64,
    inner_radius: f64,
    frames: Option<&mut Vec<RelaxFrame>>,
) {
    let ranges = molecule_ranges(particles);
    if ranges.len() < 2 {
        return;
    }
    let (outer, inner) = ranges.split_at(outer_lipids.min(ranges.len()));
    let shells: Vec<(&[std::ops::Range<usize>], f64)> =
        [(outer, outer_radius), (inner, inner_radius)]
            .into_iter()
            .filter(|(s, _)| !s.is_empty())
            .collect();
    let radii: Vec<f64> = particles
        .iter()
        .map(|p| bead_packing_radius(force_field, &p.bead))
        .collect();
    let mut positions: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();

    // When recording, project each lipid (its first bead, about the origin) to 2D
    // with a Lambert azimuthal equal-area map, so a uniform spread over the sphere
    // reads as uniform. Outer shell to the upper panel, inner to the lower.
    let mut record;
    let observer: Option<pack::StepObserver> = if let Some(frames) = frames {
        let outer = outer.to_vec();
        let inner = inner.to_vec();
        let shell_points = |pos: &[[f64; 3]], ranges: &[std::ops::Range<usize>]| {
            ranges
                .iter()
                .map(|r| lambert_azimuthal(pos[r.start]))
                .collect()
        };
        record = move |pos: &[[f64; 3]]| {
            frames.push(RelaxFrame {
                upper: shell_points(pos, &outer),
                lower: shell_points(pos, &inner),
            });
        };
        Some(&mut record)
    } else {
        None
    };

    pack::equalize_on_sphere(
        &mut positions,
        &radii,
        &shells,
        [0.0, 0.0, 0.0],
        DEFAULT_PACK_ITERATIONS,
        DEFAULT_PACK_TOLERANCE,
        observer,
    );
    for (p, pos) in particles.iter_mut().zip(positions) {
        p.position = pos;
    }
}

/// Project a point (relative to the sphere centre) to 2D with a Lambert azimuthal
/// equal-area map about the +z axis. The whole sphere lands in the disc of radius
/// 2, and equal areas on the sphere stay equal in the plane.
fn lambert_azimuthal(p: [f64; 3]) -> [f64; 2] {
    let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-9);
    let (x, y, z) = (p[0] / r, p[1] / r, p[2] / r);
    let k = (2.0 / (1.0 + z).max(1e-6)).sqrt();
    [k * x, k * y]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::membrane::molecule_ranges;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn pure(lipid: &str) -> Composition {
        Composition::from_weights([(lipid.to_string(), 1.0)]).unwrap()
    }

    fn vesicle(radius: VesicleRadius) -> Membrane {
        let ff = ForceField::isolf();
        let leaflet = pure("POPC");
        let mut rng = StdRng::seed_from_u64(1);
        build_vesicle(
            &ff,
            "vesicle",
            &leaflet,
            &leaflet,
            radius,
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap()
    }

    #[test]
    fn outer_shell_holds_more_lipids_than_inner() {
        let membrane = vesicle(VesicleRadius::Inner(8.0));
        // Split molecules by the radius of their head bead relative to the centre.
        let centre = 0.5 * membrane.box_size[0];
        let (mut inner, mut outer) = (0, 0);
        for range in molecule_ranges(&membrane.particles) {
            let head = membrane.particles[range.start].position;
            let r = ((head[0] - centre).powi(2)
                + (head[1] - centre).powi(2)
                + (head[2] - centre).powi(2))
            .sqrt();
            // Heads sit on either the inner (~8 nm) or outer (~8 + thickness) shell.
            if r < 8.0 + 2.0 {
                inner += 1
            } else {
                outer += 1
            }
        }
        assert!(outer > inner);
        assert_eq!(inner + outer, membrane.total_lipids());
    }

    #[test]
    fn inner_and_outer_radius_give_a_consistent_geometry() {
        // Building from the inner radius, then from the resulting outer radius,
        // yields the same lipid counts.
        let from_inner = vesicle(VesicleRadius::Inner(10.0));
        let outer_nm = 0.5 * from_inner.box_size[0] - crate::membrane::DEFAULT_PADDING;
        let from_outer = vesicle(VesicleRadius::Outer(outer_nm));
        assert_eq!(from_inner.total_lipids(), from_outer.total_lipids());
    }

    #[test]
    fn outer_radius_smaller_than_thickness_is_rejected() {
        let ff = ForceField::isolf();
        let leaflet = pure("POPC");
        let mut rng = StdRng::seed_from_u64(0);
        let error = build_vesicle(
            &ff,
            "v",
            &leaflet,
            &leaflet,
            VesicleRadius::Outer(2.0),
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap_err();
        assert!(matches!(error, BuildError::VesicleTooSmall { .. }));
    }

    #[test]
    fn beads_sit_between_the_two_shells() {
        let membrane = vesicle(VesicleRadius::Inner(10.0));
        let centre = 0.5 * membrane.box_size[0];
        let radii: Vec<f64> = membrane
            .particles
            .iter()
            .map(|p| {
                ((p.position[0] - centre).powi(2)
                    + (p.position[1] - centre).powi(2)
                    + (p.position[2] - centre).powi(2))
                .sqrt()
            })
            .collect();
        let min = radii.iter().cloned().fold(f64::MAX, f64::min);
        let max = radii.iter().cloned().fold(f64::MIN, f64::max);
        // Every bead lies in the bilayer shell, between the lumen and the outer
        // surface; nothing collapses to the centre.
        assert!(min > 5.0);
        assert!(max < 10.0 + 6.0);
    }

    /// A small vesicle, for the fast relaxation checks below.
    fn small_vesicle(seed: u64) -> Membrane {
        let ff = ForceField::isolf();
        let leaflet = pure("POPC");
        let mut rng = StdRng::seed_from_u64(seed);
        build_vesicle(
            &ff,
            "vesicle",
            &leaflet,
            &leaflet,
            VesicleRadius::Inner(4.0),
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap()
    }

    #[test]
    fn placement_is_reproducible_for_a_fixed_seed() {
        // Same seed ⇒ identical random placement ⇒ identical relaxed vesicle.
        assert_eq!(small_vesicle(7).particles, small_vesicle(7).particles);
    }

    #[test]
    fn relaxation_leaves_no_bead_overlap() {
        // After the spherical CVT + de-clash, no two beads of different lipids sit
        // closer than a bead contact distance.
        let v = small_vesicle(7);
        let p = &v.particles;
        let mut min_r = f64::INFINITY;
        for a in 0..p.len() {
            for b in (a + 1)..p.len() {
                if p[a].molecule_id == p[b].molecule_id {
                    continue;
                }
                let d2: f64 = p[a]
                    .position
                    .iter()
                    .zip(&p[b].position)
                    .map(|(x, y)| (x - y) * (x - y))
                    .sum();
                min_r = min_r.min(d2.sqrt());
            }
        }
        assert!(
            min_r > 0.4,
            "min bead distance {min_r} nm indicates overlap"
        );
    }

    #[test]
    fn both_radii_build_with_a_sufficient_gap() {
        let ff = ForceField::isolf();
        let leaflet = pure("POPC");
        let mut rng = StdRng::seed_from_u64(1);
        // inner 3, outer 9: a 6 nm gap comfortably holds a POPC bilayer.
        let v = build_vesicle(
            &ff,
            "v",
            &leaflet,
            &leaflet,
            VesicleRadius::Both {
                inner: 3.0,
                outer: 9.0,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap();
        assert!(v.total_lipids() > 0);
        assert!(v.box_size[0] >= 2.0 * 9.0, "box spans the outer radius");
    }

    #[test]
    fn both_radii_too_thin_is_rejected() {
        let ff = ForceField::isolf();
        let leaflet = pure("POPC");
        let mut rng = StdRng::seed_from_u64(0);
        // inner 10, outer 11: a 1 nm gap is far too thin for a POPC bilayer.
        let err = build_vesicle(
            &ff,
            "v",
            &leaflet,
            &leaflet,
            VesicleRadius::Both {
                inner: 10.0,
                outer: 11.0,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap_err();
        assert!(matches!(err, BuildError::VesicleTooThin { .. }));
    }
}
