//! Construction of a coarse-grained lipid tubule (a cylindrical bilayer,
//! periodic along `z`).
//!
//! Lipid heads are placed over an outer and an inner cylindrical shell (uniformly
//! or phase-separated into axial rings or circular patches per
//! [`BuildOptions::phase`]), then relaxed to a uniform, clash-free layout by a
//! per-shell cylindrical centroidal-Voronoi (Lloyd) pass and cylindrical de-clash
//! ([`crate::pack::equalize_on_cylinder`]). Each lipid's beads are stacked
//! **radially** away from or toward the cylinder axis (`z`), so the two tail
//! regions meet near the mid-radius.
//!
//! The resulting box spans `[2·R_outer + 2·padding, 2·R_outer + 2·padding, length]`:
//! periodic along `z` (an infinite membrane tube) with vacuum padding in `x` and
//! `y` well past the force-field cutoff.

use std::f64::consts::TAU;

use rand::Rng;

use crate::composition::Composition;
use crate::error::BuildError;
use crate::force_field::{ForceField, POTENTIAL_MINIMUM_FACTOR};
use crate::membrane::{
    BuildOptions, DEFAULT_PACK_ITERATIONS, DEFAULT_PACK_TOLERANCE, ENERGY_PACKING_AREA_FACTOR,
    Membrane, Particle, PlaneDomainLayout, RelaxFrame, aggregate_counts, canonicalize, clash_floor,
    domain_area_fractions,
};
use crate::pack;
use crate::vesicle::max_chain_length;
use crate::voronoi::wrap;

/// The tubule radius the user fixed. With one radius the other is derived from
/// the bilayer thickness; with both, the radial gap must hold the bilayer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TubuleRadius {
    /// Inner (lumen) radius in nm.
    Inner(f64),
    /// Outer radius in nm.
    Outer(f64),
    /// Both the inner (lumen) and outer radii in nm.
    Both { inner: f64, outer: f64 },
}

/// Build a lipid tubule from outer and inner shell compositions, a radius spec,
/// and a periodic axial length (nm).
///
/// `rng` seeds the random placement; pass a seeded generator for reproducible
/// output.
///
/// # Errors
///
/// Returns a [`BuildError`] if a composition names an unknown lipid or bead,
/// `length` is not positive and finite ([`InvalidLength`](BuildError::InvalidLength)),
/// or the radii cannot hold the bilayer
/// ([`TubuleTooSmall`](BuildError::TubuleTooSmall) or
/// [`TubuleTooThin`](BuildError::TubuleTooThin)).
#[allow(clippy::too_many_arguments)]
pub fn build_tubule(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: TubuleRadius,
    length: f64,
    options: &BuildOptions,
    rng: &mut impl Rng,
) -> Result<Membrane, BuildError> {
    build_system(
        force_field,
        name,
        outer,
        inner,
        radius,
        length,
        options,
        &[],
        rng,
        None,
    )
    .map(crate::embed::EmbeddedSystem::into_membrane)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Shell {
    Outer,
    Inner,
}

/// Compute `(lipid_count, cvt_radius)` for a cylindrical shell of head radius
/// `head_radius` and axial length `length`.
///
/// Each bead `k` sits at radius `r_k` with minimal-energy footprint
/// `A_k = 1.08 · (2^(1/6) · σ_k)²` on a cylinder of area `2π · r_k · length`.
/// Sizing the shell by the tightest bead (`max_k A_k / r_k`) prevents outer-tail
/// compression on tightly curved tubules while reducing to `2π · R · L / spacing²`
/// as `R → ∞`. Running cylindrical CVT at that bottleneck radius keeps the
/// densest bead layer isotropically spaced.
fn cylinder_shell_plan(
    force_field: &ForceField,
    composition: &Composition,
    head_radius: f64,
    length: f64,
    shell: Shell,
) -> Result<(usize, f64), BuildError> {
    let mut max_demand = 0.0_f64;
    let mut bottleneck_radius = head_radius;
    for component in composition.components() {
        let lipid = force_field
            .lipid(&component.name)
            .ok_or_else(|| BuildError::UnknownLipid(component.name.clone()))?;
        let mut r = head_radius;
        let last = lipid.beads.len() - 1;
        for (k, bead_name) in lipid.beads.iter().enumerate() {
            let bt = force_field
                .bead_type(bead_name)
                .ok_or_else(|| BuildError::UnknownBead(bead_name.clone()))?;
            let d = POTENTIAL_MINIMUM_FACTOR * bt.interaction.sigma();
            let demand = ENERGY_PACKING_AREA_FACTOR * d * d / r;
            if demand > max_demand {
                max_demand = demand;
                bottleneck_radius = r;
            }
            if k < last {
                r += match shell {
                    Shell::Outer => -lipid.bonds[k].length,
                    Shell::Inner => lipid.bonds[k].length,
                };
            }
        }
    }
    if max_demand <= 0.0 {
        return Err(BuildError::EmptyComposition);
    }
    let count = (TAU * length / max_demand).round().max(1.0) as usize;
    Ok((count, bottleneck_radius))
}

struct CylinderBeadTemplate<'a> {
    name: &'a str,
    charge: f64,
    mass: f64,
    packing_radius: f64,
    radius: f64,
}

#[allow(clippy::too_many_arguments)]
fn place_cylinder_shell(
    force_field: &ForceField,
    composition: &Composition,
    radius: f64,
    length: f64,
    count: usize,
    shell: Shell,
    options: &BuildOptions,
    rng: &mut impl Rng,
    buf: &mut crate::membrane::LipidBuffers,
) -> Result<Vec<usize>, BuildError> {
    let counts = composition.partition(count);
    if count == 0 {
        return Ok(counts);
    }

    let domain_fracs = domain_area_fractions(composition, &counts, force_field);
    // Unrolled coordinates: axis 0 = circumference s ∈ [0, 2πr), axis 1 = z ∈ [0, length).
    // Stripes run as axial rings along z (stripe_axis = 1); a raft centres at θ=0, z=L/2.
    let layout = PlaneDomainLayout::new(
        options.phase,
        &domain_fracs,
        [TAU * radius, length],
        1,
        [0.0, 0.5 * length],
        options.anti_register && shell == Shell::Inner,
    );

    let mut templates: Vec<(&str, Vec<CylinderBeadTemplate<'_>>, usize, usize)> =
        Vec::with_capacity(counts.len());
    let mut extra_beads = 0usize;
    for (idx, (component, &n)) in composition.components().iter().zip(&counts).enumerate() {
        if n == 0 {
            continue;
        }
        let lipid = force_field
            .lipid(&component.name)
            .ok_or_else(|| BuildError::UnknownLipid(component.name.clone()))?;
        let mut r = radius;
        let last = lipid.beads.len() - 1;
        let mut beads = Vec::with_capacity(lipid.beads.len());
        for (k, bead_name) in lipid.beads.iter().enumerate() {
            let bt = force_field
                .bead_type(bead_name)
                .ok_or_else(|| BuildError::UnknownBead(bead_name.clone()))?;
            beads.push(CylinderBeadTemplate {
                name: bt.name.as_str(),
                charge: bt.charge,
                mass: bt.mass,
                packing_radius: POTENTIAL_MINIMUM_FACTOR * bt.interaction.sigma() * 0.5,
                radius: r,
            });
            if k < last {
                r += match shell {
                    Shell::Outer => -lipid.bonds[k].length,
                    Shell::Inner => lipid.bonds[k].length,
                };
            }
        }
        extra_beads += beads.len() * n;
        templates.push((lipid.name.as_str(), beads, n, composition.domain_of(idx)));
    }

    buf.particles.reserve(extra_beads);
    buf.radii.reserve(extra_beads);
    buf.ranges.reserve(count);

    let inv_r = 1.0 / radius;
    for (lipid_name, beads, n, domain) in &templates {
        for _ in 0..*n {
            buf.molecule_id += 1;
            let mol_start = buf.particles.len();
            let [s, z] = layout.sample(*domain, rng);
            let (sin_t, cos_t) = (s * inv_r).sin_cos();
            for b in beads {
                buf.atom_id += 1;
                buf.particles.push(Particle {
                    residue: (*lipid_name).to_string(),
                    bead: b.name.to_string(),
                    molecule_id: buf.molecule_id,
                    atom_id: buf.atom_id,
                    position: [cos_t * b.radius, sin_t * b.radius, z],
                    charge: b.charge,
                    mass: b.mass,
                });
                buf.radii.push(b.packing_radius);
            }
            buf.ranges.push(mol_start..buf.particles.len());
        }
    }

    Ok(counts)
}

/// Unified cylindrical-tubule builder supporting zero, one, or many embedded proteins.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_system(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: TubuleRadius,
    length: f64,
    options: &BuildOptions,
    proteins: &[crate::embed::ProteinSpec],
    rng: &mut impl Rng,
    frames: Option<&mut Vec<RelaxFrame>>,
) -> Result<crate::embed::EmbeddedSystem, BuildError> {
    if !length.is_finite() || length <= 0.0 {
        return Err(BuildError::InvalidLength(length));
    }

    let gap = clash_floor(force_field, outer)?.max(clash_floor(force_field, inner)?);
    let thickness =
        max_chain_length(force_field, outer)? + max_chain_length(force_field, inner)? + gap;

    let (inner_radius, outer_radius) = match radius {
        TubuleRadius::Inner(r) => {
            if !r.is_finite() || r <= 0.0 {
                return Err(BuildError::TubuleTooSmall {
                    outer: r,
                    thickness,
                });
            }
            (r, r + thickness)
        }
        TubuleRadius::Outer(r) => {
            if !r.is_finite() || r <= thickness {
                return Err(BuildError::TubuleTooSmall {
                    outer: r,
                    thickness,
                });
            }
            (r - thickness, r)
        }
        TubuleRadius::Both { inner, outer } => {
            if !inner.is_finite() || inner <= 0.0 {
                return Err(BuildError::TubuleTooSmall {
                    outer: inner,
                    thickness,
                });
            }
            if !outer.is_finite() || outer - inner < thickness {
                return Err(BuildError::TubuleTooThin {
                    inner,
                    outer,
                    minimum: thickness,
                });
            }
            (inner, outer)
        }
    };

    let r_mid = 0.5 * (inner_radius + outer_radius);
    let placed = if proteins.is_empty() {
        None
    } else {
        Some(crate::embed::place_proteins_tubule(
            proteins, r_mid, length,
        )?)
    };

    // Size each shell from its tightest bead layer so curved packing stays
    // clash-free at every radius and approaches area / spacing² as R → ∞.
    let (outer_count, outer_cvt_radius) =
        cylinder_shell_plan(force_field, outer, outer_radius, length, Shell::Outer)?;
    let (inner_count, inner_cvt_radius) =
        cylinder_shell_plan(force_field, inner, inner_radius, length, Shell::Inner)?;

    let mut buf = crate::membrane::LipidBuffers::default();
    let outer_counts = place_cylinder_shell(
        force_field,
        outer,
        outer_radius,
        length,
        outer_count,
        Shell::Outer,
        options,
        rng,
        &mut buf,
    )?;
    let mut outer_lipids = buf.molecule_id;
    let inner_counts = place_cylinder_shell(
        force_field,
        inner,
        inner_radius,
        length,
        inner_count,
        Shell::Inner,
        options,
        rng,
        &mut buf,
    )?;

    let metric = crate::embed::CarveMetric::Cylinder { length };
    let carved = match &placed {
        Some(p) => crate::embed::carve_lipids(
            &mut buf.particles,
            &mut buf.radii,
            &mut buf.ranges,
            &mut outer_lipids,
            &p.obstacles,
            0.0,
            metric,
        ),
        None => 0,
    };

    let obstacles: &[pack::Obstacle] = placed.as_ref().map_or(&[], |p| &p.obstacles);
    relax_tubule(
        &mut buf,
        outer_lipids,
        obstacles,
        [outer_radius, outer_cvt_radius, inner_cvt_radius],
        length,
        frames,
    );

    let pruned = match &placed {
        Some(p) => crate::embed::prune_clashing_lipids(
            &mut buf.particles,
            &mut buf.radii,
            &mut buf.ranges,
            &mut outer_lipids,
            &p.raw_obstacles,
            metric,
        ),
        None => 0,
    };

    let lipid_counts = if carved + pruned > 0 {
        crate::embed::count_remaining_lipids(&buf.particles, outer, inner)
    } else {
        aggregate_counts(outer, &outer_counts, inner, &inner_counts)
    };
    let mut particles = canonicalize(buf.particles, &lipid_counts);

    let mut max_r = outer_radius;
    if let Some(p) = &placed {
        for pos in p.iter_positions() {
            max_r = max_r.max(pos[0].hypot(pos[1]));
        }
    }
    // Centre the cylinder axis in the padded x/y cross-section; z stays in [0, length).
    let box_xy = 2.0 * max_r + 2.0 * options.padding;
    let half_xy = 0.5 * box_xy;
    for particle in &mut particles {
        particle.position[0] += half_xy;
        particle.position[1] += half_xy;
    }
    if let Some(mut placed) = placed {
        placed.translate([half_xy, half_xy, 0.0]);
        Ok(crate::embed::assemble_system(
            name,
            placed,
            particles,
            lipid_counts,
            [box_xy, box_xy, length],
            carved,
            pruned,
        ))
    } else {
        Ok(crate::embed::EmbeddedSystem::from_membrane(Membrane {
            name: name.to_string(),
            box_size: [box_xy, box_xy, length],
            lipid_counts,
            particles,
        }))
    }
}

fn relax_tubule(
    buf: &mut crate::membrane::LipidBuffers,
    outer_lipids: usize,
    obstacles: &[pack::Obstacle],
    radii_spec: [f64; 3],
    length: f64,
    frames: Option<&mut Vec<RelaxFrame>>,
) {
    let [outer_radius, outer_cvt_radius, inner_cvt_radius] = radii_spec;
    let particles = &mut buf.particles;
    let radii = &buf.radii;
    let ranges = &buf.ranges;
    if ranges.len() < 2 {
        return;
    }
    let (outer, inner) = ranges.split_at(outer_lipids.min(ranges.len()));
    let shells: Vec<(&[std::ops::Range<usize>], f64)> =
        [(outer, outer_cvt_radius), (inner, inner_cvt_radius)]
            .into_iter()
            .filter(|(s, _)| !s.is_empty())
            .collect();
    let mut positions: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
    match frames {
        Some(out) => {
            out.push(project_cylinder_heads(
                &positions,
                ranges,
                outer_lipids,
                outer_radius,
                length,
            ));
            let mut record = |pos: &[[f64; 3]]| {
                out.push(project_cylinder_heads(
                    pos,
                    ranges,
                    outer_lipids,
                    outer_radius,
                    length,
                ));
            };
            pack::equalize_on_cylinder_with_obstacles(
                &mut positions,
                radii,
                &shells,
                obstacles,
                length,
                DEFAULT_PACK_ITERATIONS,
                DEFAULT_PACK_TOLERANCE,
                Some(&mut record),
            );
        }
        None => {
            pack::equalize_on_cylinder_with_obstacles(
                &mut positions,
                radii,
                &shells,
                obstacles,
                length,
                DEFAULT_PACK_ITERATIONS,
                DEFAULT_PACK_TOLERANCE,
                None,
            );
        }
    }
    for (p, pos) in particles.iter_mut().zip(positions) {
        p.position = pos;
    }
}

/// Project every lipid's head bead onto the unrolled cylinder `[0, 2π·R_outer) × [0, length)`.
fn project_cylinder_heads(
    positions: &[[f64; 3]],
    ranges: &[std::ops::Range<usize>],
    outer_lipids: usize,
    outer_radius: f64,
    length: f64,
) -> RelaxFrame {
    let project = |slice: &[std::ops::Range<usize>]| {
        slice
            .iter()
            .map(|r| {
                let p = positions[r.start];
                let theta = p[1].atan2(p[0]).rem_euclid(TAU);
                [theta * outer_radius, wrap(p[2], length)]
            })
            .collect()
    };
    let (outer, inner) = ranges.split_at(outer_lipids.min(ranges.len()));
    RelaxFrame {
        upper: project(outer),
        lower: project(inner),
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    use super::*;
    use crate::membrane::Phase;

    fn popc() -> Composition {
        Composition::parse("POPC").unwrap()
    }

    #[test]
    fn tubule_outer_shell_holds_more_lipids_than_inner_and_matches_box() {
        let ff = ForceField::isolf();
        let mut rng = StdRng::seed_from_u64(1);
        let opts = BuildOptions::default();
        let t = build_tubule(
            &ff,
            "t",
            &popc(),
            &popc(),
            TubuleRadius::Outer(10.0),
            20.0,
            &opts,
            &mut rng,
        )
        .unwrap();

        let expected_xy = 2.0 * 10.0 + 2.0 * opts.padding;
        assert!((t.box_size[0] - expected_xy).abs() < 1e-9);
        assert!((t.box_size[1] - expected_xy).abs() < 1e-9);
        assert!((t.box_size[2] - 20.0).abs() < 1e-9);

        let mid = 0.5 * expected_xy;
        let head_radii: Vec<f64> = t
            .particles
            .chunk_by(|a, b| a.molecule_id == b.molecule_id)
            .map(|mol| (mol[0].position[0] - mid).hypot(mol[0].position[1] - mid))
            .collect();
        let r_max = head_radii.iter().copied().fold(0.0_f64, f64::max);
        let r_min = head_radii.iter().copied().fold(f64::MAX, f64::min);
        assert!((r_max - 10.0).abs() < 1e-6);
        assert!(r_min < r_max - 2.0);

        let cut = 0.5 * (r_min + r_max);
        let outer_n = head_radii.iter().filter(|&&r| r > cut).count();
        let inner_n = head_radii.iter().filter(|&&r| r <= cut).count();
        assert!(outer_n > inner_n, "outer {outer_n} <= inner {inner_n}");
    }

    #[test]
    fn tubule_rejects_invalid_geometry() {
        let ff = ForceField::isolf();
        let mut rng = StdRng::seed_from_u64(1);
        let opts = BuildOptions::default();

        assert!(matches!(
            build_tubule(
                &ff,
                "t",
                &popc(),
                &popc(),
                TubuleRadius::Outer(2.0),
                15.0,
                &opts,
                &mut rng
            ),
            Err(BuildError::TubuleTooSmall { .. })
        ));
        assert!(matches!(
            build_tubule(
                &ff,
                "t",
                &popc(),
                &popc(),
                TubuleRadius::Both {
                    inner: 8.0,
                    outer: 9.0
                },
                15.0,
                &opts,
                &mut rng
            ),
            Err(BuildError::TubuleTooThin { .. })
        ));
        assert!(matches!(
            build_tubule(
                &ff,
                "t",
                &popc(),
                &popc(),
                TubuleRadius::Outer(10.0),
                0.0,
                &opts,
                &mut rng
            ),
            Err(BuildError::InvalidLength(_))
        ));
    }

    #[test]
    fn tubule_has_no_inter_lipid_clashes_across_periodic_z() {
        let ff = ForceField::isolf();
        let mut rng = StdRng::seed_from_u64(42);
        let length = 15.0;
        let t = build_tubule(
            &ff,
            "t",
            &popc(),
            &popc(),
            TubuleRadius::Outer(8.0),
            length,
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap();

        let half_z = 0.5 * length;
        let radii: Vec<f64> = t
            .particles
            .iter()
            .map(|p| {
                let bt = ff.bead_type(&p.bead).unwrap();
                POTENTIAL_MINIMUM_FACTOR * bt.interaction.sigma() * 0.5
            })
            .collect();
        let mut max_overlap = 0.0_f64;
        for i in 0..t.particles.len() {
            for j in (i + 1)..t.particles.len() {
                if t.particles[i].molecule_id == t.particles[j].molecule_id {
                    continue;
                }
                let dx = t.particles[i].position[0] - t.particles[j].position[0];
                let dy = t.particles[i].position[1] - t.particles[j].position[1];
                let mut dz = t.particles[i].position[2] - t.particles[j].position[2];
                if dz > half_z {
                    dz -= length;
                } else if dz < -half_z {
                    dz += length;
                }
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                max_overlap = max_overlap.max(radii[i] + radii[j] - d);
            }
        }
        assert!(max_overlap < 0.05, "residual tubule overlap {max_overlap}");
    }

    #[test]
    fn phase_stripe_and_raft_separate_domains_on_tubule() {
        let ff = ForceField::isolf();
        let comp = Composition::parse("DPPC=1,POPC=1").unwrap();

        // Axial rings along z.
        let mut rng = StdRng::seed_from_u64(9);
        let striped = build_tubule(
            &ff,
            "t",
            &comp,
            &comp,
            TubuleRadius::Outer(9.0),
            24.0,
            &BuildOptions {
                phase: Phase::Stripe,
                ..BuildOptions::default()
            },
            &mut rng,
        )
        .unwrap();
        let mean_z = |res: &str| {
            let zs: Vec<f64> = striped
                .particles
                .chunk_by(|a, b| a.molecule_id == b.molecule_id)
                .filter(|m| m[0].residue == res)
                .map(|m| m[0].position[2])
                .collect();
            zs.iter().sum::<f64>() / zs.len() as f64
        };
        assert!(
            (mean_z("DPPC") - mean_z("POPC")).abs() > 8.0,
            "axial stripes did not separate: DPPC z={}, POPC z={}",
            mean_z("DPPC"),
            mean_z("POPC")
        );

        // Circular raft on +x side vs anti-registered inner shell.
        let mut rng = StdRng::seed_from_u64(9);
        let raft_comp = Composition::parse("DPPC=1,POPC=3").unwrap();
        let raft = build_tubule(
            &ff,
            "t",
            &raft_comp,
            &raft_comp,
            TubuleRadius::Outer(10.0),
            24.0,
            &BuildOptions {
                phase: Phase::Raft,
                anti_register: true,
                ..BuildOptions::default()
            },
            &mut rng,
        )
        .unwrap();
        let mid_xy = 0.5 * raft.box_size[0];
        let mut outer_dppc_x = 0.0;
        let mut outer_dppc_n = 0usize;
        let mut inner_dppc_x = 0.0;
        let mut inner_dppc_n = 0usize;
        let r_cut = 10.0 - 2.0;
        for mol in raft
            .particles
            .chunk_by(|a, b| a.molecule_id == b.molecule_id)
            .filter(|m| m[0].residue == "DPPC")
        {
            let dx = mol[0].position[0] - mid_xy;
            let dy = mol[0].position[1] - mid_xy;
            if dx.hypot(dy) > r_cut {
                outer_dppc_x += dx;
                outer_dppc_n += 1;
            } else {
                inner_dppc_x += dx;
                inner_dppc_n += 1;
            }
        }
        let outer_mean = outer_dppc_x / outer_dppc_n as f64;
        let inner_mean = inner_dppc_x / inner_dppc_n as f64;
        assert!(
            outer_mean > 3.0 && inner_mean < -2.0,
            "anti-registered tubule rafts should sit on opposite sides: outer={outer_mean}, inner={inner_mean}"
        );
    }
}
