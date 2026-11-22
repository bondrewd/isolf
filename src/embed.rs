//! Embedding coarse-grained proteins into membranes, vesicles, and tubules.
//!
//! Proteins keep their input PDB/mmCIF orientation by default (OPM convention:
//! bilayer midplane at `z = 0`, outward normal `+z`), with optional overrides
//! for vertical placement (`z=center` or offset), lateral coordinates, tilt,
//! axial rotation, copy count, and lipid clearance.
//!
//! For each geometry, embedding follows a four-stage pipeline:
//! 1. **Place** protein copies on the bilayer mid-surface (flat `z = 0`,
//!    spherical `R_mid`, or cylindrical `R_mid`) and construct excluded-volume
//!    [`Obstacle`] spheres (`r_prot + clearance`).
//! 2. **Carve** initial lipids whose beads overlap any obstacle.
//! 3. **Relax** remaining lipids with obstacle-aware pinned-site Voronoi CVT
//!    and soft-sphere de-clash ([`crate::pack`]).
//! 4. **Prune** any lipid still trapped in a residual clash against a protein
//!    bead, then assemble a unified [`EmbeddedSystem`] with matching `.gro`,
//!    `.top`, `.psf`, `.pdb`, `.crd`, `.cif`, `.itp`, and `.vmd` outputs.

use std::f64::consts::TAU;
use std::fmt::Write;
use std::ops::Range;

use rand::Rng;

use crate::composition::Composition;
use crate::error::BuildError;
use crate::force_field::ForceField;
use crate::membrane::{BuildOptions, Membrane, Particle, RelaxFrame, Sizing, molecule_ranges};
use crate::pack::Obstacle;
use crate::protein::{Protein, protein_segname};
use crate::tubule::TubuleRadius;
use crate::vesicle::VesicleRadius;
use crate::voronoi::wrap;

/// Golden angle in radians (`π * (3 - √5)`), used for uniform Fibonacci-sphere
/// and helical multi-copy layouts.
const GOLDEN_ANGLE: f64 = 2.399_963_229_728_653;

/// Default extra lipid-carving clearance around protein Cα beads (nm).
pub const DEFAULT_CLEARANCE: f64 = 0.15;

/// Residual overlap threshold (nm) above which a trapped lipid is pruned after
/// de-clash.
const POST_RELAX_PRUNE_TOL: f64 = 0.05;

/// How a protein's vertical (`z`, or radial on curved bilayers) position is set.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ZPlacement {
    /// Preserve the input structure's `z` coordinates relative to the bilayer
    /// midplane `z = 0` (OPM convention, default).
    #[default]
    Keep,
    /// Shift the protein's centroid `z` to the bilayer midplane `z = 0`.
    Center,
}

/// Placement options for embedding one protein species into a lipid geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinPlacement {
    /// Vertical/radial baseline (`Keep` preserves OPM `z`, `Center` centers centroid).
    pub z_mode: ZPlacement,
    /// Additional shift along the membrane normal (nm, default `0.0`).
    pub z_offset: f64,
    /// Optional lateral x shift (nm; circumferential shift on tubule; tangent on vesicle).
    pub x: Option<f64>,
    /// Optional lateral y shift (nm; axial shift on tubule; tangent on vesicle).
    pub y: Option<f64>,
    /// Tilt away from the membrane normal (degrees, default `0.0`).
    pub tilt_deg: f64,
    /// Rotation about the membrane normal (degrees, default `0.0`).
    pub rot_deg: f64,
    /// Number of copies to embed (default `1`).
    pub copies: usize,
    /// Extra steric clearance added to each protein bead's radius during lipid
    /// carving and de-clash (nm, default [`DEFAULT_CLEARANCE`]).
    pub clearance: f64,
}

impl Default for ProteinPlacement {
    fn default() -> Self {
        Self {
            z_mode: ZPlacement::Keep,
            z_offset: 0.0,
            x: None,
            y: None,
            tilt_deg: 0.0,
            rot_deg: 0.0,
            copies: 1,
            clearance: DEFAULT_CLEARANCE,
        }
    }
}

/// A coarse-grained [`Protein`] paired with its [`ProteinPlacement`] instructions.
#[derive(Debug, Clone)]
pub struct ProteinSpec {
    pub protein: Protein,
    pub placement: ProteinPlacement,
}

/// A built lipid + protein system (supports zero, one, or many embedded proteins).
#[derive(Debug, Clone, PartialEq)]
pub struct EmbeddedSystem {
    /// System title written into output headers.
    pub name: String,
    /// Unique protein templates (one `.itp` file written per template).
    pub protein_templates: Vec<Protein>,
    /// Instance counts per template in `[ molecules ]` order.
    pub protein_counts: Vec<(String, usize)>,
    /// All placed protein copies in `[ molecules ]` order (matching `particles`).
    pub proteins: Vec<Protein>,
    /// Number of lipids carved before relaxation to accommodate proteins.
    pub carved_lipids: usize,
    /// Number of trapped clashing lipids pruned after relaxation.
    pub pruned_lipids: usize,
    /// Non-fatal warnings generated during protein placement/embedding.
    pub warnings: Vec<String>,
    /// Combined particle list: all protein beads first, then canonicalized lipid beads.
    pub particles: Vec<Particle>,
    /// Periodic box dimensions (nm).
    pub box_size: [f64; 3],
    /// Remaining lipid counts per species in canonical order.
    pub lipid_counts: Vec<(String, usize)>,
    /// Starting index of lipid particles in `particles` (`0` when no proteins).
    pub lipid_particle_start: usize,
}

impl EmbeddedSystem {
    /// Wrap a lipid-only [`Membrane`] with zero proteins.
    pub fn from_membrane(membrane: Membrane) -> Self {
        Self {
            name: membrane.name,
            protein_templates: Vec::new(),
            protein_counts: Vec::new(),
            proteins: Vec::new(),
            carved_lipids: 0,
            pruned_lipids: 0,
            warnings: Vec::new(),
            particles: membrane.particles,
            box_size: membrane.box_size,
            lipid_counts: membrane.lipid_counts,
            lipid_particle_start: 0,
        }
    }

    /// True if at least one protein is embedded.
    pub fn has_proteins(&self) -> bool {
        !self.proteins.is_empty()
    }

    /// Total number of remaining lipid molecules.
    pub fn total_lipids(&self) -> usize {
        self.lipid_counts.iter().map(|(_, c)| c).sum()
    }

    /// Total number of embedded protein instances.
    pub fn total_proteins(&self) -> usize {
        self.proteins.len()
    }

    /// Total number of molecule entries in `[ molecules ]` (proteins + lipids).
    pub fn total_molecules(&self) -> usize {
        self.total_proteins() + self.total_lipids()
    }

    /// Slice containing only the lipid particles.
    pub fn lipid_particles(&self) -> &[Particle] {
        &self.particles[self.lipid_particle_start..]
    }

    /// Translate every particle (and stored protein copy) by `shift` (nm).
    pub fn translate(&mut self, shift: [f64; 3]) {
        for p in &mut self.particles {
            p.position[0] += shift[0];
            p.position[1] += shift[1];
            p.position[2] += shift[2];
        }
        for prot in &mut self.proteins {
            prot.translate(shift);
        }
    }

    /// Centre the system inside `[0, box_size]`.
    ///
    /// For lipid-only systems this matches [`Membrane::center_in_box`]; when
    /// proteins are present, the bilayer/system centre is shifted cleanly by
    /// `box_size / 2` from its origin-centred build frame (preserving OPM `z = 0`
    /// when shifted back by `-box_size / 2` for `--center origin`).
    pub fn center_in_box(&mut self) {
        if self.particles.is_empty() {
            return;
        }
        if !self.has_proteins() {
            let n = self.particles.len() as f64;
            let mut sum = [0.0; 3];
            for p in &self.particles {
                sum[0] += p.position[0];
                sum[1] += p.position[1];
                sum[2] += p.position[2];
            }
            self.translate([
                self.box_size[0] / 2.0 - sum[0] / n,
                self.box_size[1] / 2.0 - sum[1] / n,
                self.box_size[2] / 2.0 - sum[2] / n,
            ]);
            return;
        }

        // Shift using the lipid particles' centroid so asymmetric extracellular
        // protein domains do not drag the bilayer midplane off-center.
        let lipids = self.lipid_particles();
        let src = if lipids.is_empty() {
            self.particles.as_slice()
        } else {
            lipids
        };
        let n = src.len() as f64;
        let mut sum = [0.0; 3];
        for p in src {
            sum[0] += p.position[0];
            sum[1] += p.position[1];
            sum[2] += p.position[2];
        }
        let mut shift = [
            self.box_size[0] / 2.0 - sum[0] / n,
            self.box_size[1] / 2.0 - sum[1] / n,
            self.box_size[2] / 2.0 - sum[2] / n,
        ];
        // Ensure no protruding protein bead lands outside [0, box_z].
        let (mut z_lo, mut z_hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for p in &self.particles {
            let z = p.position[2] + shift[2];
            z_lo = z_lo.min(z);
            z_hi = z_hi.max(z);
        }
        if z_lo < 0.0 || z_hi > self.box_size[2] {
            let mid = 0.5 * (z_lo + z_hi);
            shift[2] += 0.5 * self.box_size[2] - mid;
        }
        self.translate(shift);
    }

    /// Render the combined system as a GROMACS `.gro` coordinate file.
    pub fn to_gro(&self) -> String {
        let mut s = String::with_capacity(self.particles.len() * 72 + 128);
        let _ = writeln!(s, "{}, t = {:>16.3} ", self.name, 0.0);
        let _ = writeln!(s, "{:>12} ", self.particles.len());
        for p in &self.particles {
            let _ = writeln!(
                s,
                "{:>5}{:>5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}  0.0000  0.0000  0.0000 ",
                p.molecule_id % 100000,
                p.residue,
                p.bead,
                p.atom_id % 100000,
                p.position[0],
                p.position[1],
                p.position[2],
            );
        }
        let _ = writeln!(
            s,
            "{:>15.4}{:>15.4}{:>15.4} ",
            self.box_size[0], self.box_size[1], self.box_size[2]
        );
        s.push('\n');
        s
    }

    /// Render the combined system as a GENESIS `.top` topology file.
    pub fn to_top(&self) -> String {
        let mut s = String::new();
        if self.has_proteins() {
            let _ = writeln!(s, "; AICG2+ protein force field parameters");
            let _ = writeln!(s, "#include \"./aicg2p.itp\"");
        }
        let _ = writeln!(s, "; Implicit Solvent Lipid Forcefield (iSoLF)");
        let _ = writeln!(s, "#include \"./isolf.itp\"");
        if self.has_proteins() {
            let _ = writeln!(s, "; Protein topologies");
            for tmpl in &self.protein_templates {
                let _ = writeln!(s, "#include \"./{}.itp\"", tmpl.name);
            }
        }
        let _ = writeln!(s);
        let _ = writeln!(s, "[ system ]");
        let _ = writeln!(s, "{}", self.name);
        let _ = writeln!(s);
        let _ = writeln!(s, "[ molecules ]");
        for (name, count) in &self.protein_counts {
            let _ = writeln!(s, "{name} {count}");
        }
        for (lipid, count) in &self.lipid_counts {
            let _ = writeln!(s, "{lipid} {count}");
        }
        let _ = writeln!(s);
        let _ = writeln!(s, "[ cg_ele_chain_pairs ]");
        let total = self.total_molecules();
        let _ = writeln!(s, "ON 1 - {total} : 1 - {total}");
        s
    }

    /// Render the combined system as a chain-aware CHARMM `.psf` file.
    pub fn to_psf(&self) -> String {
        if !self.has_proteins() {
            return crate::psf::render(&self.particles, &self.name);
        }
        let mut chains: Vec<(String, Range<usize>)> = Vec::new();
        let mut offset = 0usize;
        let mut prot_chain_counter = 0usize;
        for prot in &self.proteins {
            for ch in &prot.chains {
                let seg = protein_segname(prot_chain_counter, &ch.seg_name);
                let start = offset + ch.range.start;
                let end = offset + ch.range.end;
                chains.push((seg, start..end));
                prot_chain_counter += 1;
            }
            offset += prot.beads.len();
        }
        let lipid_ranges = molecule_ranges(self.lipid_particles());
        for r in lipid_ranges {
            chains.push(("MEMB".to_string(), (offset + r.start)..(offset + r.end)));
        }
        crate::psf::render_chains(&self.particles, &chains, &self.name)
    }

    /// Render the combined system as a `.pdb` coordinate file (coordinates in Å).
    pub fn to_pdb(&self) -> String {
        if !self.has_proteins() {
            return crate::pdb::render(&self.particles, &self.name);
        }
        let mut s = String::with_capacity(self.particles.len() * 82 + 256);
        let _ = writeln!(s, "TITLE     {:<70}", self.name);
        let mut offset = 0usize;
        let mut chain_counter = 0usize;
        for prot in &self.proteins {
            for ch in &prot.chains {
                let cid = if chain_counter < 26 {
                    (b'A' + chain_counter as u8) as char
                } else {
                    'Z'
                };
                for idx in ch.range.clone() {
                    let p = &self.particles[offset + idx];
                    let _ = writeln!(
                        s,
                        "ATOM  {:>5} {:>4} {:<4}{}{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00               ",
                        p.atom_id % 100000,
                        p.bead,
                        p.residue,
                        cid,
                        p.molecule_id % 10000,
                        p.position[0] * 10.0,
                        p.position[1] * 10.0,
                        p.position[2] * 10.0,
                    );
                }
                s.push_str("TER\n");
                chain_counter += 1;
            }
            offset += prot.beads.len();
        }
        for p in self.lipid_particles() {
            let _ = writeln!(
                s,
                "ATOM  {:>5} {:>4} {:<4}M{:>4}    {:>8.3}{:>8.3}{:>8.3}  0.00  0.00               ",
                p.atom_id % 100000,
                p.bead,
                p.residue,
                p.molecule_id % 10000,
                p.position[0] * 10.0,
                p.position[1] * 10.0,
                p.position[2] * 10.0,
            );
        }
        s.push_str("TER\nEND\n");
        s
    }

    /// Render the combined system as a CHARMM `.crd` coordinate file.
    pub fn to_crd(&self) -> String {
        let num_res = self.lipid_particle_start + self.total_lipids();
        crate::crd::render(&self.particles, &self.name, num_res)
    }

    /// Render the combined system as an mmCIF `.cif` coordinate file.
    pub fn to_cif(&self) -> String {
        crate::cif::render(&self.particles, &self.name)
    }
}

// ==================== Public geometry entry points ====================

/// Build a flat lipid bilayer with embedded proteins.
///
/// # Errors
///
/// Returns [`BuildError`] on invalid lipid compositions, sizing, or protein
/// placement options.
#[allow(clippy::too_many_arguments)]
pub fn embed_membrane(
    force_field: &ForceField,
    name: &str,
    upper: &Composition,
    lower: &Composition,
    sizing: Sizing,
    options: &BuildOptions,
    proteins: &[ProteinSpec],
    rng: &mut impl Rng,
) -> Result<EmbeddedSystem, BuildError> {
    crate::membrane::build_with_proteins(
        force_field,
        name,
        upper,
        lower,
        sizing,
        options,
        proteins,
        rng,
        None,
    )
}

/// Like [`embed_membrane`], also recording relaxation frames for `--gif`.
///
/// # Errors
///
/// The same as [`embed_membrane`].
#[allow(clippy::too_many_arguments)]
pub fn embed_membrane_recorded(
    force_field: &ForceField,
    name: &str,
    upper: &Composition,
    lower: &Composition,
    sizing: Sizing,
    options: &BuildOptions,
    proteins: &[ProteinSpec],
    rng: &mut impl Rng,
) -> Result<(EmbeddedSystem, Vec<RelaxFrame>), BuildError> {
    let mut frames = Vec::new();
    let sys = crate::membrane::build_with_proteins(
        force_field,
        name,
        upper,
        lower,
        sizing,
        options,
        proteins,
        rng,
        Some(&mut frames),
    )?;
    Ok((sys, frames))
}

/// Build a spherical vesicle with embedded proteins.
///
/// # Errors
///
/// Returns [`BuildError`] on invalid lipid compositions, vesicle radii, or
/// protein placement options.
#[allow(clippy::too_many_arguments)]
pub fn embed_vesicle(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: VesicleRadius,
    options: &BuildOptions,
    proteins: &[ProteinSpec],
    rng: &mut impl Rng,
) -> Result<EmbeddedSystem, BuildError> {
    crate::vesicle::build_with_proteins(
        force_field,
        name,
        outer,
        inner,
        radius,
        options,
        proteins,
        rng,
        None,
    )
}

/// Like [`embed_vesicle`], also recording relaxation frames for `--gif`.
///
/// # Errors
///
/// The same as [`embed_vesicle`].
#[allow(clippy::too_many_arguments)]
pub fn embed_vesicle_recorded(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: VesicleRadius,
    options: &BuildOptions,
    proteins: &[ProteinSpec],
    rng: &mut impl Rng,
) -> Result<(EmbeddedSystem, Vec<RelaxFrame>), BuildError> {
    let mut frames = Vec::new();
    let sys = crate::vesicle::build_with_proteins(
        force_field,
        name,
        outer,
        inner,
        radius,
        options,
        proteins,
        rng,
        Some(&mut frames),
    )?;
    Ok((sys, frames))
}

/// Build a cylindrical tubule with embedded proteins.
///
/// # Errors
///
/// Returns [`BuildError`] on invalid lipid compositions, tubule dimensions, or
/// protein placement options.
#[allow(clippy::too_many_arguments)]
pub fn embed_tubule(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: TubuleRadius,
    length: f64,
    options: &BuildOptions,
    proteins: &[ProteinSpec],
    rng: &mut impl Rng,
) -> Result<EmbeddedSystem, BuildError> {
    crate::tubule::build_with_proteins(
        force_field,
        name,
        outer,
        inner,
        radius,
        length,
        options,
        proteins,
        rng,
        None,
    )
}

/// Like [`embed_tubule`], also recording relaxation frames for `--gif`.
///
/// # Errors
///
/// The same as [`embed_tubule`].
#[allow(clippy::too_many_arguments)]
pub fn embed_tubule_recorded(
    force_field: &ForceField,
    name: &str,
    outer: &Composition,
    inner: &Composition,
    radius: TubuleRadius,
    length: f64,
    options: &BuildOptions,
    proteins: &[ProteinSpec],
    rng: &mut impl Rng,
) -> Result<(EmbeddedSystem, Vec<RelaxFrame>), BuildError> {
    let mut frames = Vec::new();
    let sys = crate::tubule::build_with_proteins(
        force_field,
        name,
        outer,
        inner,
        radius,
        length,
        options,
        proteins,
        rng,
        Some(&mut frames),
    )?;
    Ok((sys, frames))
}

// ==================== Internal placement, carving & assembly ====================

/// Validate all placement parameters up front.
fn validate_specs(specs: &[ProteinSpec]) -> Result<(), BuildError> {
    for spec in specs {
        let p = &spec.placement;
        if p.copies == 0 {
            return Err(BuildError::InvalidPlacement(format!(
                "protein '{}': copies must be at least 1",
                spec.protein.name
            )));
        }
        if !p.clearance.is_finite() || p.clearance < 0.0 {
            return Err(BuildError::InvalidPlacement(format!(
                "protein '{}': clearance must be non-negative (got {})",
                spec.protein.name, p.clearance
            )));
        }
        if !p.z_offset.is_finite() || !p.tilt_deg.is_finite() || !p.rot_deg.is_finite() {
            return Err(BuildError::InvalidPlacement(format!(
                "protein '{}': placement offsets and angles must be finite",
                spec.protein.name
            )));
        }
        if let Some(x) = p.x
            && !x.is_finite()
        {
            return Err(BuildError::InvalidPlacement(format!(
                "protein '{}': x must be finite",
                spec.protein.name
            )));
        }
        if let Some(y) = p.y
            && !y.is_finite()
        {
            return Err(BuildError::InvalidPlacement(format!(
                "protein '{}': y must be finite",
                spec.protein.name
            )));
        }
    }
    Ok(())
}

/// Deduplicate protein templates and group placed copies in canonical `[ molecules ]`
/// order, returning `(templates, counts, instances_ordered)`.
pub(crate) struct PlacedProteins {
    pub templates: Vec<Protein>,
    pub counts: Vec<(String, usize)>,
    pub instances: Vec<Protein>,
    pub obstacles: Vec<Obstacle>,
    pub raw_obstacles: Vec<Obstacle>,
    pub warnings: Vec<String>,
}

/// Compare two protein topologies (ignoring Cartesian coordinates) so identical
/// templates can share one `.itp` and one `[ molecules ]` entry.
fn same_topology(a: &Protein, b: &Protein) -> bool {
    a.name == b.name
        && a.model == b.model
        && a.chains == b.chains
        && a.bonds == b.bonds
        && a.angles == b.angles
        && a.dihedrals == b.dihedrals
        && a.pairs == b.pairs
        && a.beads.len() == b.beads.len()
        && a.beads.iter().zip(&b.beads).all(|(x, y)| {
            x.residue == y.residue
                && x.atom_type == y.atom_type
                && (x.charge - y.charge).abs() < 1e-9
                && (x.mass - y.mass).abs() < 1e-9
        })
}

/// Center protein `(x, y)` on `(0, 0)`, apply `z_mode` + `z_offset`, then tilt
/// and rotate around `(0, 0, z_offset)`.
fn orient_local(template: &Protein, p: &ProteinPlacement, extra_rot_deg: f64) -> Protein {
    let mut prot = template.clone();
    let c = prot.centroid();
    let dz = match p.z_mode {
        ZPlacement::Keep => p.z_offset,
        ZPlacement::Center => -c[2] + p.z_offset,
    };
    prot.translate([-c[0], -c[1], dz]);
    let pivot = [0.0, 0.0, p.z_offset];
    if p.tilt_deg.abs() > 1e-9 {
        prot.rotate_axis(pivot, [0.0, 1.0, 0.0], p.tilt_deg.to_radians());
    }
    let rot = p.rot_deg + extra_rot_deg;
    if rot.abs() > 1e-9 {
        prot.rotate_axis(pivot, [0.0, 0.0, 1.0], rot.to_radians());
    }
    prot
}

/// Check if an oriented local protein sits far from the bilayer midplane `z = 0`.
fn check_midplane_warning(
    prot: &Protein,
    placement: &ProteinPlacement,
    warnings: &mut Vec<String>,
) {
    if placement.z_mode != ZPlacement::Keep {
        return;
    }
    let (lo, hi) = prot.bounding_box();
    if hi[2] < -4.0 || lo[2] > 4.0 {
        warnings.push(format!(
            "protein '{}' sits at z = [{:.1}, {:.1}] nm, far from the bilayer midplane z = 0 \
             (pass z=center to center it in the membrane)",
            prot.name, lo[2], hi[2]
        ));
    }
}

/// Group placed copies into deduplicated templates with unique `[ moleculetype ]` names.
fn finalize_placed(
    specs: &[ProteinSpec],
    per_spec_copies: Vec<Vec<Protein>>,
    warnings: Vec<String>,
) -> PlacedProteins {
    let mut templates: Vec<Protein> = Vec::new();
    let mut grouped: Vec<Vec<(Protein, f64)>> = Vec::new();
    let mut used_names: Vec<String> = Vec::new();

    for (spec, copies) in specs.iter().zip(per_spec_copies) {
        let clearance = spec.placement.clearance;
        if let Some(pos) = templates
            .iter()
            .position(|t| same_topology(t, &spec.protein))
        {
            let tmpl_name = templates[pos].name.clone();
            for mut inst in copies {
                inst.name = tmpl_name.clone();
                grouped[pos].push((inst, clearance));
            }
        } else {
            let mut tmpl = spec.protein.clone();
            if used_names.contains(&tmpl.name) {
                let base = tmpl.name.clone();
                let mut suffix = 2usize;
                while used_names.contains(&format!("{base}_{suffix}")) {
                    suffix += 1;
                }
                tmpl.name = format!("{base}_{suffix}");
            }
            let tmpl_name = tmpl.name.clone();
            used_names.push(tmpl_name.clone());
            templates.push(tmpl);
            let mut list = Vec::with_capacity(copies.len());
            for mut inst in copies {
                inst.name = tmpl_name.clone();
                list.push((inst, clearance));
            }
            grouped.push(list);
        }
    }

    let mut counts = Vec::with_capacity(templates.len());
    let mut instances = Vec::new();
    let mut obstacles = Vec::new();
    let mut raw_obstacles = Vec::new();
    for (tmpl, list) in templates.iter().zip(grouped) {
        counts.push((tmpl.name.clone(), list.len()));
        for (inst, clearance) in list {
            for b in &inst.beads {
                obstacles.push(Obstacle {
                    pos: b.position,
                    radius: b.radius + clearance,
                });
                raw_obstacles.push(Obstacle {
                    pos: b.position,
                    radius: b.radius,
                });
            }
            instances.push(inst);
        }
    }

    PlacedProteins {
        templates,
        counts,
        instances,
        obstacles,
        raw_obstacles,
        warnings,
    }
}

/// Place protein copies in a flat membrane box `[0, box_x) × [0, box_y)` with
/// bilayer midplane at `z = 0`.
pub(crate) fn place_proteins_membrane(
    specs: &[ProteinSpec],
    box_x: f64,
    box_y: f64,
) -> Result<PlacedProteins, BuildError> {
    validate_specs(specs)?;
    let mut warnings = Vec::new();

    // Count how many instances need automatic lateral grid placement.
    let auto_total: usize = specs
        .iter()
        .map(|s| {
            if s.placement.x.is_some() || s.placement.y.is_some() {
                if s.placement.copies > 1 {
                    s.placement.copies
                } else {
                    0
                }
            } else {
                s.placement.copies
            }
        })
        .sum();
    let cols = (auto_total as f64 * (box_x / box_y.max(1e-6)))
        .sqrt()
        .ceil()
        .max(1.0) as usize;
    let rows = auto_total.div_ceil(cols).max(1);
    let mut auto_idx = 0usize;

    let mut per_spec_copies = Vec::with_capacity(specs.len());
    for spec in specs {
        let mut copies = Vec::with_capacity(spec.placement.copies);
        let base_local = orient_local(&spec.protein, &spec.placement, 0.0);
        check_midplane_warning(&base_local, &spec.placement, &mut warnings);

        let (lo, hi) = base_local.bounding_box();
        if (hi[0] - lo[0]) > box_x || (hi[1] - lo[1]) > box_y {
            warnings.push(format!(
                "protein '{}' lateral span ({:.1} x {:.1} nm) exceeds membrane box ({:.1} x {:.1} nm)",
                spec.protein.name,
                hi[0] - lo[0],
                hi[1] - lo[1],
                box_x,
                box_y
            ));
        }

        for k in 0..spec.placement.copies {
            let extra_rot = if spec.placement.copies > 1 {
                (k as f64) * 360.0 / (spec.placement.copies as f64)
            } else {
                0.0
            };
            let mut inst = if k == 0 && extra_rot == 0.0 {
                base_local.clone()
            } else {
                orient_local(&spec.protein, &spec.placement, extra_rot)
            };

            let (cx, cy) = if spec.placement.copies == 1
                && (spec.placement.x.is_some() || spec.placement.y.is_some())
            {
                (
                    wrap(0.5 * box_x + spec.placement.x.unwrap_or(0.0), box_x),
                    wrap(0.5 * box_y + spec.placement.y.unwrap_or(0.0), box_y),
                )
            } else if auto_total == 1 {
                (0.5 * box_x, 0.5 * box_y)
            } else {
                let col = auto_idx % cols;
                let row = auto_idx / cols;
                auto_idx += 1;
                let gx = (col as f64 + 0.5) / cols as f64 * box_x;
                let gy = (row as f64 + 0.5) / rows as f64 * box_y;
                (
                    wrap(gx + spec.placement.x.unwrap_or(0.0), box_x),
                    wrap(gy + spec.placement.y.unwrap_or(0.0), box_y),
                )
            };
            inst.translate([cx, cy, 0.0]);
            copies.push(inst);
        }
        per_spec_copies.push(copies);
    }

    Ok(finalize_placed(specs, per_spec_copies, warnings))
}

/// Place protein copies radially on a spherical vesicle mid-surface of radius `r_mid`.
pub(crate) fn place_proteins_vesicle(
    specs: &[ProteinSpec],
    r_mid: f64,
) -> Result<PlacedProteins, BuildError> {
    validate_specs(specs)?;
    let mut warnings = Vec::new();
    let total_copies: usize = specs.iter().map(|s| s.placement.copies).sum();
    let mut global_k = 0usize;

    let mut per_spec_copies = Vec::with_capacity(specs.len());
    for spec in specs {
        let mut copies = Vec::with_capacity(spec.placement.copies);
        let base_local = orient_local(&spec.protein, &spec.placement, 0.0);
        check_midplane_warning(&base_local, &spec.placement, &mut warnings);

        for k in 0..spec.placement.copies {
            let mut inst = if k == 0 {
                base_local.clone()
            } else {
                orient_local(&spec.protein, &spec.placement, 0.0)
            };

            let mut normal = if total_copies == 1 {
                [0.0, 0.0, 1.0]
            } else {
                let z = 1.0 - (2.0 * global_k as f64 + 1.0) / total_copies as f64;
                let r = (1.0 - z * z).max(0.0).sqrt();
                let phi = global_k as f64 * GOLDEN_ANGLE;
                let (sp, cp) = phi.sin_cos();
                [r * cp, r * sp, z]
            };
            global_k += 1;

            // Optional tangent displacement (x, y) in nm along the sphere surface.
            let (sx, sy) = (
                spec.placement.x.unwrap_or(0.0),
                spec.placement.y.unwrap_or(0.0),
            );
            let arc = sx.hypot(sy);
            if arc > 1e-9 && r_mid > 1e-6 {
                let (u0, v0, _) = crate::protein::orthonormal_frame(normal);
                let dir = [
                    (sx * u0[0] + sy * v0[0]) / arc,
                    (sx * u0[1] + sy * v0[1]) / arc,
                    (sx * u0[2] + sy * v0[2]) / arc,
                ];
                let (sa, ca) = (arc / r_mid).sin_cos();
                normal = [
                    ca * normal[0] + sa * dir[0],
                    ca * normal[1] + sa * dir[1],
                    ca * normal[2] + sa * dir[2],
                ];
            }

            let (u, v, n) = crate::protein::orthonormal_frame(normal);
            let origin = [r_mid * n[0], r_mid * n[1], r_mid * n[2]];
            inst.apply_frame(origin, u, v, n);
            copies.push(inst);
        }
        per_spec_copies.push(copies);
    }

    Ok(finalize_placed(specs, per_spec_copies, warnings))
}

/// Place protein copies radially on a cylindrical tubule mid-surface of radius
/// `r_mid` and periodic axial length `length`.
pub(crate) fn place_proteins_tubule(
    specs: &[ProteinSpec],
    r_mid: f64,
    length: f64,
) -> Result<PlacedProteins, BuildError> {
    validate_specs(specs)?;
    let mut warnings = Vec::new();
    let total_copies: usize = specs.iter().map(|s| s.placement.copies).sum();
    let mut global_k = 0usize;

    let mut per_spec_copies = Vec::with_capacity(specs.len());
    for spec in specs {
        let mut copies = Vec::with_capacity(spec.placement.copies);
        let base_local = orient_local(&spec.protein, &spec.placement, 0.0);
        check_midplane_warning(&base_local, &spec.placement, &mut warnings);

        for k in 0..spec.placement.copies {
            let mut inst = if k == 0 {
                base_local.clone()
            } else {
                orient_local(&spec.protein, &spec.placement, 0.0)
            };

            let (base_theta, base_z) = if total_copies == 1 {
                (0.0, 0.5 * length)
            } else {
                let z = (global_k as f64 + 0.5) / total_copies as f64 * length;
                let th = (global_k as f64 * GOLDEN_ANGLE).rem_euclid(TAU);
                (th, z)
            };
            global_k += 1;

            let theta = base_theta + spec.placement.x.unwrap_or(0.0) / r_mid.max(1e-6);
            let z_cyl = wrap(base_z + spec.placement.y.unwrap_or(0.0), length);
            let (st, ct) = theta.sin_cos();
            let n = [ct, st, 0.0];
            let u = [-st, ct, 0.0];
            let v = [0.0, 0.0, 1.0];
            let origin = [r_mid * ct, r_mid * st, z_cyl];
            inst.apply_frame(origin, u, v, n);
            copies.push(inst);
        }
        per_spec_copies.push(copies);
    }

    Ok(finalize_placed(specs, per_spec_copies, warnings))
}

/// Remove every lipid molecule in `ranges` that has any bead overlapping any
/// obstacle according to `dist_sq(p_lipid, p_obs) < (r_lipid + r_obs - margin)²`.
/// Compacts `particles`, `radii`, `ranges`, and updates `first_leaflet_lipids`.
/// Returns the number of removed lipids.
pub(crate) fn carve_lipids(
    particles: &mut Vec<Particle>,
    radii: &mut Vec<f64>,
    ranges: &mut Vec<Range<usize>>,
    first_leaflet_lipids: &mut usize,
    obstacles: &[Obstacle],
    margin: f64,
    dist_sq: impl Fn([f64; 3], [f64; 3]) -> f64,
) -> usize {
    if obstacles.is_empty() || ranges.is_empty() {
        return 0;
    }

    let mut keep_lipid = vec![true; ranges.len()];
    let mut removed = 0usize;
    for (li, range) in ranges.iter().enumerate() {
        'lipid: for k in range.clone() {
            let p = particles[k].position;
            let rl = radii[k];
            for o in obstacles {
                let cut = (rl + o.radius - margin).max(0.0);
                if cut > 0.0 && dist_sq(p, o.pos) < cut * cut {
                    keep_lipid[li] = false;
                    removed += 1;
                    break 'lipid;
                }
            }
        }
    }
    if removed == 0 {
        return 0;
    }

    let old_first = *first_leaflet_lipids;
    let mut new_particles = Vec::with_capacity(particles.len());
    let mut new_radii = Vec::with_capacity(radii.len());
    let mut new_ranges = Vec::with_capacity(ranges.len() - removed);
    let mut new_first = 0usize;
    for (li, range) in ranges.iter().enumerate() {
        if !keep_lipid[li] {
            continue;
        }
        if li < old_first {
            new_first += 1;
        }
        let start = new_particles.len();
        new_particles.extend(particles[range.clone()].iter().cloned());
        new_radii.extend_from_slice(&radii[range.clone()]);
        new_ranges.push(start..new_particles.len());
    }
    *particles = new_particles;
    *radii = new_radii;
    *ranges = new_ranges;
    *first_leaflet_lipids = new_first;
    removed
}

/// Remove any trapped lipids still clashing with raw protein beads after de-clash.
pub(crate) fn prune_clashing_lipids(
    particles: &mut Vec<Particle>,
    radii: &mut Vec<f64>,
    ranges: &mut Vec<Range<usize>>,
    first_leaflet_lipids: &mut usize,
    raw_obstacles: &[Obstacle],
    dist_sq: impl Fn([f64; 3], [f64; 3]) -> f64,
) -> usize {
    carve_lipids(
        particles,
        radii,
        ranges,
        first_leaflet_lipids,
        raw_obstacles,
        POST_RELAX_PRUNE_TOL,
        dist_sq,
    )
}

/// Count remaining lipid molecules per species in canonical order (matching
/// `upper`/`lower` first-appearance order, excluding zero-count species).
pub(crate) fn count_remaining_lipids(
    particles: &[Particle],
    upper: &Composition,
    lower: &Composition,
) -> Vec<(String, usize)> {
    let mut order: Vec<String> = Vec::new();
    for comp in [upper, lower] {
        for c in comp.components() {
            if !order.contains(&c.name) {
                order.push(c.name.clone());
            }
        }
    }
    let mut counts = vec![0usize; order.len()];
    for mol in particles.chunk_by(|a, b| a.molecule_id == b.molecule_id) {
        if let Some(pos) = order.iter().position(|s| *s == mol[0].residue) {
            counts[pos] += 1;
        }
    }
    order
        .into_iter()
        .zip(counts)
        .filter(|(_, c)| *c > 0)
        .collect()
}

/// Combine placed proteins and canonicalized lipids into a single
/// [`EmbeddedSystem`] with continuous 1-based residue/molecule and atom IDs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_system(
    name: &str,
    placed: PlacedProteins,
    mut lipid_particles: Vec<Particle>,
    lipid_counts: Vec<(String, usize)>,
    box_size: [f64; 3],
    carved_lipids: usize,
    pruned_lipids: usize,
) -> EmbeddedSystem {
    let mut particles = Vec::new();
    let mut next_res_id = 1usize;
    let mut next_atom_id = 1usize;

    for prot in &placed.instances {
        let prot_parts = prot.to_particles(next_res_id, next_atom_id);
        next_res_id += prot.beads.len();
        next_atom_id += prot.beads.len();
        particles.extend(prot_parts);
    }

    let lipid_particle_start = particles.len();
    let res_offset = next_res_id - 1;
    let atom_offset = next_atom_id - 1;
    for p in &mut lipid_particles {
        p.molecule_id += res_offset;
        p.atom_id += atom_offset;
    }
    particles.extend(lipid_particles);

    EmbeddedSystem {
        name: name.to_string(),
        protein_templates: placed.templates,
        protein_counts: placed.counts,
        proteins: placed.instances,
        carved_lipids,
        pruned_lipids,
        warnings: placed.warnings,
        particles,
        box_size,
        lipid_counts,
        lipid_particle_start,
    }
}

/// Periodic 2D (x/y) + Euclidean z squared distance for flat bilayers.
#[inline]
pub(crate) fn dist_sq_plane(a: [f64; 3], b: [f64; 3], box_x: f64, box_y: f64) -> f64 {
    let mut dx = (a[0] - b[0]).abs();
    if dx > 0.5 * box_x {
        dx = box_x - dx;
    }
    let mut dy = (a[1] - b[1]).abs();
    if dy > 0.5 * box_y {
        dy = box_y - dy;
    }
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Euclidean 3D squared distance for vesicles.
#[inline]
pub(crate) fn dist_sq_3d(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Periodic z + Euclidean x/y squared distance for tubules.
#[inline]
pub(crate) fn dist_sq_cylinder(a: [f64; 3], b: [f64; 3], length: f64) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let mut dz = (a[2] - b[2]).abs();
    if dz > 0.5 * length {
        dz = length - dz;
    }
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    use super::*;
    use crate::protein::aicg2p::ProteinOptions;
    use crate::protein::structure::parse_pdb;

    /// A 22-residue transmembrane α-helix aligned along z, centered at z = 0
    /// (spans z ≈ -1.6 to +1.6 nm, matching an OPM-oriented TM helix).
    fn tm_helix_protein() -> Protein {
        let mut pdb = String::new();
        let mut serial = 1;
        for i in 0..22 {
            let idx = i + 1;
            let theta = i as f64 * (100.0_f64).to_radians();
            let r = 2.3;
            let z = (i as f64 - 10.5) * 1.5; // -15.75 to +15.75 Å
            let ca = [r * theta.cos(), r * theta.sin(), z];
            let n = [r * (theta - 0.3).cos(), r * (theta - 0.3).sin(), z - 0.6];
            let c = [r * (theta + 0.3).cos(), r * (theta + 0.3).sin(), z + 0.6];
            let o = [
                (r + 0.4) * (theta + 0.3).cos(),
                (r + 0.4) * (theta + 0.3).sin(),
                z + 1.3,
            ];
            for (aname, p) in [("N", n), ("CA", ca), ("C", c), ("O", o)] {
                let _ = writeln!(
                    pdb,
                    "ATOM  {:>5} {:<4} LEU A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00",
                    serial, aname, idx, p[0], p[1], p[2]
                );
                serial += 1;
            }
        }
        pdb.push_str("TER\nEND\n");
        let structure = parse_pdb(&pdb).unwrap();
        Protein::from_structure("tmhelix_cg", &structure, &ProteinOptions::default()).unwrap()
    }

    #[test]
    fn embeds_tm_helix_in_membrane_without_protein_lipid_clashes() {
        let ff = ForceField::isolf();
        let comp = Composition::parse("POPC").unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        let spec = ProteinSpec {
            protein: tm_helix_protein(),
            placement: ProteinPlacement::default(),
        };
        let sys = embed_membrane(
            &ff,
            "membrane",
            &comp,
            &comp,
            Sizing::Count {
                upper: 64,
                lower: 64,
            },
            &BuildOptions::default(),
            &[spec],
            &mut rng,
        )
        .unwrap();

        assert!(sys.has_proteins());
        assert_eq!(sys.total_proteins(), 1);
        assert_eq!(sys.lipid_particle_start, 22);
        assert!(
            sys.carved_lipids > 0,
            "expected carved lipids around TM helix"
        );
        assert!(sys.total_lipids() < 128);
        assert!(
            sys.warnings.is_empty(),
            "unexpected warnings: {:?}",
            sys.warnings
        );

        // Verify zero steric clash between protein and lipid beads.
        let (bx, by) = (sys.box_size[0], sys.box_size[1]);
        let mut min_dist = f64::INFINITY;
        for pb in &sys.proteins[0].beads {
            for lp in sys.lipid_particles() {
                let d = dist_sq_plane(pb.position, lp.position, bx, by).sqrt();
                min_dist = min_dist.min(d);
            }
        }
        assert!(
            min_dist > 0.45,
            "protein-lipid bead clash detected: min dist = {min_dist:.3} nm"
        );

        // Verify combined topology & PSF consistency.
        let top = sys.to_top();
        assert!(top.contains("#include \"./aicg2p.itp\""));
        assert!(top.contains("#include \"./isolf.itp\""));
        assert!(top.contains("#include \"./tmhelix_cg.itp\""));
        assert!(top.contains("tmhelix_cg 1"));
        assert!(top.contains(&format!(
            "ON 1 - {} : 1 - {}",
            sys.total_molecules(),
            sys.total_molecules()
        )));

        let psf = sys.to_psf();
        assert!(psf.contains("PROA"));
        assert!(psf.contains("MEMB"));
        assert!(psf.contains(&format!("{:>8} !NATOM", sys.particles.len())));
    }

    #[test]
    fn embeds_multiple_copies_in_vesicle_and_tubule() {
        let ff = ForceField::isolf();
        let comp = Composition::parse("POPC").unwrap();
        let spec = ProteinSpec {
            protein: tm_helix_protein(),
            placement: ProteinPlacement {
                copies: 3,
                ..ProteinPlacement::default()
            },
        };

        let mut rng = StdRng::seed_from_u64(7);
        let ves = embed_vesicle(
            &ff,
            "vesicle",
            &comp,
            &comp,
            VesicleRadius::Outer(8.0),
            &BuildOptions::default(),
            std::slice::from_ref(&spec),
            &mut rng,
        )
        .unwrap();
        assert_eq!(ves.total_proteins(), 3);
        assert_eq!(ves.protein_counts, vec![("tmhelix_cg".to_string(), 3)]);
        assert!(ves.carved_lipids > 0);
        for prot in &ves.proteins {
            for pb in &prot.beads {
                for lp in ves.lipid_particles() {
                    let d = dist_sq_3d(pb.position, lp.position).sqrt();
                    assert!(d > 0.45, "vesicle protein-lipid clash: d={d:.3}");
                }
            }
        }

        let mut rng = StdRng::seed_from_u64(7);
        let tub = embed_tubule(
            &ff,
            "tubule",
            &comp,
            &comp,
            TubuleRadius::Outer(8.0),
            16.0,
            &BuildOptions::default(),
            &[spec],
            &mut rng,
        )
        .unwrap();
        assert_eq!(tub.total_proteins(), 3);
        assert!(tub.carved_lipids > 0);
        for prot in &tub.proteins {
            for pb in &prot.beads {
                for lp in tub.lipid_particles() {
                    let d = dist_sq_cylinder(pb.position, lp.position, 16.0).sqrt();
                    assert!(d > 0.45, "tubule protein-lipid clash: d={d:.3}");
                }
            }
        }
    }
}
