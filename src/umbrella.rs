//! 1D bilayer umbrella-sampling setup (independent windows and REUS),
//! GENESIS `RESTR_CVS` CV extraction, and MBAR/BO campaign utilities.
//!
//! Supports both single amino-acid probes (`aa:TRP`, `all` — for iSoLFv2
//! lipid–protein parameterization) and multi-residue peptide/protein probes
//! from PDB/mmCIF files (for insertion/binding PMFs and validation):
//!
//! - **Constant lipid count across every window** (zero lipids carved,
//!   required by MBAR/WHAM), using periodic lateral clearance search plus
//!   rigid-molecule lateral lipid de-clash for multi-bead probes.
//! - **PBC-safe bilayer centering** via mass-weighted circular mean in `z`,
//!   handling pre-equilibrated `--bilayer` `.gro` files even if the bilayer
//!   drifted across the periodic `z` boundary.
//! - **Exact atom-index probe selection** (`ano:1` or `ano:1-N`) paired with
//!   signed COM distance restraints (`function1 = DISTMASS`, `direction1 = Z2`,
//!   `U = k * (dz - z0)^2`), avoiding GROMACS per-molecule `resnr = 1`
//!   collisions and working identically for single beads and multi-residue
//!   peptides/proteins.
//! - **Manifest-free CV extraction & REUS demuxing** directly from GENESIS
//!   `atdyn` logs into `cv/{}.dat` for `mbar_analysis`.

use std::f64::consts::TAU;
use std::fmt::Write as _;
use std::path::Path;

use crate::embed::EmbeddedSystem;
use crate::error::BuildError;
use crate::force_field::{ForceField, POTENTIAL_MINIMUM_FACTOR};
use crate::membrane::{Membrane, Particle, molecule_ranges};
use crate::protein::Protein;

pub use crate::force_field::lipid_protein::{
    DEFAULT_PROBE_EPSILON, DEFAULT_PROBE_OMEGA, IsolfBeadParam, LipidProteinParams, LjPairParam,
    ProbeParamOverrides, WcaPairParam,
};

/// Default harmonic umbrella force constant `k` in kcal·mol⁻¹·Å⁻²
/// (`U = k * (dz - z0)^2`, matching GENESIS `atdyn` & `mbar_analysis`).
pub const DEFAULT_UMBRELLA_K: f64 = 0.5;

/// Default z-box length (nm) for bilayer umbrella sampling (`200 Å`, matching
/// Ugarte La Torre & Takada 2021).
pub const DEFAULT_UMBRELLA_BOX_Z: f64 = 20.0;

/// Classification of an umbrella-sampling solute probe (extensible for future
/// lipid flip-flop probes).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeKind {
    /// Single coarse-grained amino-acid bead (e.g. `TRP`, `ARG`).
    AminoAcid,
    /// Multi-bead peptide or protein coarse-grained from a PDB/mmCIF structure.
    Protein,
}

impl ProbeKind {
    /// Short machine-readable label written into `manifest.json`.
    pub fn label(self) -> &'static str {
        match self {
            Self::AminoAcid => "amino_acid",
            Self::Protein => "protein",
        }
    }
}

/// A prepared solute probe ready to place across umbrella `z` windows.
#[derive(Debug, Clone)]
pub struct UmbrellaProbe {
    /// Probe category (`AminoAcid` or `Protein`).
    pub kind: ProbeKind,
    /// Coarse-grained coordinates and topology (already oriented with any
    /// requested `tilt`/`rot` before window placement).
    pub protein: Protein,
}

impl UmbrellaProbe {
    /// Build a GENESIS `[SELECTION]` `group1` expression selecting strictly
    /// this probe's atoms (`1..=N_probe`), avoiding GROMACS `resnr = 1`
    /// collisions across lipids and working for both single beads and
    /// multi-residue peptides/proteins.
    pub fn selection(&self) -> String {
        probe_selection(self.protein.beads.len())
    }
}

/// GENESIS `[SELECTION]` expression selecting the probe atoms `1..=n_atoms`.
pub fn probe_selection(n_atoms: usize) -> String {
    if n_atoms <= 1 {
        "ano:1".to_string()
    } else {
        format!("ano:1-{n_atoms}")
    }
}

/// Parse an umbrella z-window specification (in nm):
/// - Range with step: `start..end:step` (e.g. `0.0..4.0:0.1`, inclusive)
/// - Range with count: `start..end/count` (e.g. `0.0..4.0/41`, inclusive)
/// - Explicit comma-separated list: `0.0,0.5,1.0,2.0,3.5`
///
/// # Errors
///
/// Returns [`BuildError::InvalidPlacement`] if the spec is malformed or yields
/// zero windows.
pub fn parse_z_windows(spec: &str) -> Result<Vec<f64>, BuildError> {
    let spec = spec.trim();
    if spec.is_empty() {
        return Err(BuildError::InvalidPlacement(
            "--z specification is empty".to_string(),
        ));
    }

    let parse_f = |s: &str, what: &str| -> Result<f64, BuildError> {
        s.trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| {
                BuildError::InvalidPlacement(format!("--z: invalid {what} '{}'", s.trim()))
            })
    };

    if let Some((start_str, rest)) = spec.split_once("..") {
        let start = parse_f(start_str, "start z")?;
        if let Some((end_str, step_str)) = rest.split_once(':') {
            let end = parse_f(end_str, "end z")?;
            let step = parse_f(step_str, "step")?;
            if step <= 0.0 {
                return Err(BuildError::InvalidPlacement(format!(
                    "--z: step must be positive (got {step})"
                )));
            }
            let span = end - start;
            if span < -1e-9 {
                return Err(BuildError::InvalidPlacement(format!(
                    "--z: end ({end}) must be >= start ({start})"
                )));
            }
            let n = ((span / step) + 1e-9).floor() as usize;
            let mut out = Vec::with_capacity(n + 1);
            for k in 0..=n {
                let z = clean_zero(start + (k as f64) * step);
                if z <= end + 1e-7 {
                    out.push(z);
                }
            }
            return Ok(out);
        }
        if let Some((end_str, count_str)) = rest.split_once('/') {
            let end = parse_f(end_str, "end z")?;
            let count: usize = count_str.trim().parse().map_err(|_| {
                BuildError::InvalidPlacement(format!(
                    "--z: invalid window count '{}'",
                    count_str.trim()
                ))
            })?;
            if count == 0 {
                return Err(BuildError::InvalidPlacement(
                    "--z: window count must be at least 1".to_string(),
                ));
            }
            if count == 1 {
                return Ok(vec![clean_zero(start)]);
            }
            let mut out = Vec::with_capacity(count);
            for k in 0..count {
                let frac = k as f64 / (count - 1) as f64;
                out.push(clean_zero(start + frac * (end - start)));
            }
            return Ok(out);
        }
        return Err(BuildError::InvalidPlacement(format!(
            "--z: range '{spec}' needs ':step' or '/count' (e.g. '0.0..4.0:0.1' or '0.0..4.0/41')"
        )));
    }

    let mut out = Vec::new();
    for tok in spec.split(',') {
        let tok = tok.trim();
        if !tok.is_empty() {
            out.push(clean_zero(parse_f(tok, "z coordinate")?));
        }
    }
    if out.is_empty() {
        return Err(BuildError::InvalidPlacement(
            "--z produced no windows".to_string(),
        ));
    }
    Ok(out)
}

#[inline]
fn clean_zero(v: f64) -> f64 {
    let rounded = (v * 1e6).round() / 1e6;
    if rounded == 0.0 { 0.0 } else { rounded }
}

/// Center a lipid bilayer's mass-weighted COM at `(0, 0, 0)` using PBC-safe
/// circular-mean unwrapping along `z`, and ensure `box_size[2] >= min_box_z`.
///
/// Unwrapping via circular mean first guarantees that a pre-equilibrated
/// `--bilayer` `.gro` whose bilayer drifted across the periodic `z` boundary
/// is reassembled intact at `z = 0` rather than placing the COM in vacuum.
pub fn center_bilayer_at_origin(bilayer: &mut Membrane, min_box_z: f64) {
    if bilayer.particles.is_empty() {
        return;
    }
    // 1. Unwrap z about the mass-weighted circular mean so a PBC-split bilayer
    //    becomes one contiguous slab before computing the linear COM.
    let lz = bilayer.box_size[2];
    if lz > 0.0 {
        let (mut sin_sum, mut cos_sum) = (0.0, 0.0);
        for p in &bilayer.particles {
            let m = if p.mass > 0.0 { p.mass } else { 1.0 };
            let theta = TAU * p.position[2] / lz;
            sin_sum += m * theta.sin();
            cos_sum += m * theta.cos();
        }
        if sin_sum.hypot(cos_sum) > 1e-12 {
            let z_mid = (sin_sum.atan2(cos_sum) / TAU) * lz;
            let half_z = 0.5 * lz;
            for p in &mut bilayer.particles {
                let mut dz = (p.position[2] - z_mid) % lz;
                if dz > half_z {
                    dz -= lz;
                } else if dz < -half_z {
                    dz += lz;
                }
                p.position[2] = z_mid + dz;
            }
        }
    }

    // 2. Shift exact unwrapped mass-weighted COM to (0, 0, 0).
    let mut mass_sum = 0.0;
    let mut com = [0.0; 3];
    for p in &bilayer.particles {
        let m = if p.mass > 0.0 { p.mass } else { 1.0 };
        mass_sum += m;
        com[0] += m * p.position[0];
        com[1] += m * p.position[1];
        com[2] += m * p.position[2];
    }
    let inv = 1.0 / mass_sum;
    com[0] *= inv;
    com[1] *= inv;
    com[2] *= inv;
    for p in &mut bilayer.particles {
        p.position[0] -= com[0];
        p.position[1] -= com[1];
        p.position[2] -= com[2];
    }
    if bilayer.box_size[2] < min_box_z {
        bilayer.box_size[2] = min_box_z;
    }
}

/// Build one umbrella-sampling window coordinate set with `probe`'s
/// mass-weighted COM placed at `z = z_nm` (atoms `1..=N_probe`, molecule 1)
/// and ALL lipids of `bilayer` preserved intact (molecules `2..=N_lipids+1`).
///
/// 1. Searches a 24×24 periodic lateral grid to maximize minimum bead-to-bead
///    clearance between all probe beads and lipid beads in the probe's `z` slab.
/// 2. If any overlap remains (typical when inserting a multi-bead peptide or
///    helix into the bilayer core), runs rigid-molecule 2D lateral lipid
///    de-clash so lipids part smoothly around the probe **without carving any
///    lipids** and without altering any lipid's `z` coordinate.
pub fn build_window_particles(
    bilayer: &Membrane,
    probe: &Protein,
    z_nm: f64,
    force_field: &ForceField,
) -> Vec<Particle> {
    let com = probe.center_of_mass();
    let probe_rel: Vec<([f64; 3], f64)> = probe
        .beads
        .iter()
        .map(|b| {
            (
                [
                    b.position[0] - com[0],
                    b.position[1] - com[1],
                    b.position[2] - com[2] + z_nm,
                ],
                b.radius,
            )
        })
        .collect();

    let [bx, by, _] = bilayer.box_size;
    let (px, py, best_clearance) = best_lateral_xy(bilayer, &probe_rel, bx, by, force_field);

    let n_probe = probe.beads.len();
    let mut out = Vec::with_capacity(n_probe + bilayer.particles.len());
    for (idx, (b, &(rel_pos, _))) in probe.beads.iter().zip(&probe_rel).enumerate() {
        out.push(Particle {
            residue: b.residue.clone(),
            bead: b.atom.clone(),
            molecule_id: 1,
            atom_id: idx + 1,
            position: [px + rel_pos[0], py + rel_pos[1], rel_pos[2]],
            charge: b.charge,
            mass: b.mass,
        });
    }
    for p in &bilayer.particles {
        let mut lp = p.clone();
        lp.molecule_id += 1;
        lp.atom_id += n_probe;
        out.push(lp);
    }

    if best_clearance < 0.0 {
        let probe_placed: Vec<([f64; 3], f64)> = probe_rel
            .iter()
            .map(|&(pos, r)| ([px + pos[0], py + pos[1], pos[2]], r))
            .collect();
        declash_lipids_around_probe(&mut out[n_probe..], &probe_placed, bx, by, force_field);
    }
    out
}

/// Build an [`EmbeddedSystem`] for one umbrella window so `.gro`, `.pdb`,
/// `.psf`, `.top`, and VMD scripts share one code path with `membrane -p`.
pub fn build_window_system(
    title: String,
    bilayer: &Membrane,
    probe: &Protein,
    z_nm: f64,
    force_field: &ForceField,
) -> EmbeddedSystem {
    let particles = build_window_particles(bilayer, probe, z_nm, force_field);
    EmbeddedSystem::from_probe_and_bilayer(title, probe.clone(), bilayer, particles)
}

/// Search a 24×24 periodic lateral grid over `[-bx/2, bx/2) × [-by/2, by/2)` to
/// find `(px, py)` maximizing the minimum bead-to-bead clearance between all
/// probe beads and lipid beads in the probe's vertical slab.
fn best_lateral_xy(
    bilayer: &Membrane,
    probe_rel: &[([f64; 3], f64)],
    bx: f64,
    by: f64,
    force_field: &ForceField,
) -> (f64, f64, f64) {
    const GRID: usize = 24;
    let mut z_lo = f64::INFINITY;
    let mut z_hi = f64::NEG_INFINITY;
    for &(pos, r) in probe_rel {
        z_lo = z_lo.min(pos[2] - r);
        z_hi = z_hi.max(pos[2] + r);
    }
    let slab: Vec<([f64; 3], f64)> = bilayer
        .particles
        .iter()
        .filter(|p| p.position[2] >= z_lo - 1.2 && p.position[2] <= z_hi + 1.2)
        .map(|p| {
            let r = lipid_bead_radius(&p.bead, force_field);
            (p.position, r)
        })
        .collect();
    if slab.is_empty() {
        return (0.0, 0.0, f64::INFINITY);
    }

    let half_x = 0.5 * bx;
    let half_y = 0.5 * by;
    let mut best_xy = (0.0, 0.0);
    let mut best_clearance = f64::NEG_INFINITY;
    for ix in 0..GRID {
        let px = -half_x + ((ix as f64 + 0.5) / GRID as f64) * bx;
        for iy in 0..GRID {
            let py = -half_y + ((iy as f64 + 0.5) / GRID as f64) * by;
            let mut min_clear = f64::INFINITY;
            'probe_loop: for &(ppos, pr) in probe_rel {
                let qx = px + ppos[0];
                let qy = py + ppos[1];
                let qz = ppos[2];
                for &(lpos, lr) in &slab {
                    let dz = qz - lpos[2];
                    let sum_r = pr + lr;
                    if dz.abs() - sum_r >= min_clear {
                        continue;
                    }
                    let mut dx = (qx - lpos[0]).abs() % bx;
                    if dx > half_x {
                        dx = bx - dx;
                    }
                    let mut dy = (qy - lpos[1]).abs() % by;
                    if dy > half_y {
                        dy = by - dy;
                    }
                    let d = (dx * dx + dy * dy + dz * dz).sqrt() - sum_r;
                    if d < min_clear {
                        min_clear = d;
                        if min_clear <= best_clearance {
                            break 'probe_loop;
                        }
                    }
                }
            }
            if min_clear > best_clearance {
                best_clearance = min_clear;
                best_xy = (px, py);
            }
        }
    }
    (best_xy.0, best_xy.1, best_clearance)
}

#[inline]
fn lipid_bead_radius(bead: &str, force_field: &ForceField) -> f64 {
    force_field
        .bead_type(bead)
        .map(|bt| 0.5 * POTENTIAL_MINIMUM_FACTOR * bt.interaction.sigma())
        .unwrap_or(0.35)
}

#[inline]
fn pbc_delta(d: f64, box_len: f64) -> f64 {
    let mut x = d % box_len;
    let half = 0.5 * box_len;
    if x > half {
        x -= box_len;
    } else if x < -half {
        x += box_len;
    }
    x
}

/// Rigid-molecule 2D lateral soft-sphere de-clash of lipids around a stationary
/// multi-bead probe. Preserves exact lipid count, intramolecular bonds, and all
/// `z` coordinates (so the bilayer COM stays identically at `z = 0`).
fn declash_lipids_around_probe(
    lipid_particles: &mut [Particle],
    probe: &[([f64; 3], f64)],
    bx: f64,
    by: f64,
    force_field: &ForceField,
) {
    const MAX_ITERS: usize = 60;
    const MARGIN: f64 = 0.04;
    const MAX_STEP: f64 = 0.15;

    let ranges = molecule_ranges(lipid_particles);
    let n_mols = ranges.len();
    if n_mols == 0 || probe.is_empty() {
        return;
    }
    let bead_radii: Vec<f64> = lipid_particles
        .iter()
        .map(|p| lipid_bead_radius(&p.bead, force_field))
        .collect();
    let mol_z_mean: Vec<f64> = ranges
        .iter()
        .map(|r| {
            let sum: f64 = lipid_particles[r.clone()]
                .iter()
                .map(|p| p.position[2])
                .sum();
            sum / r.len() as f64
        })
        .collect();
    let mut mol_xy: Vec<[f64; 2]> = ranges
        .iter()
        .map(|r| {
            let n = r.len() as f64;
            let (mut sx, mut sy) = (0.0, 0.0);
            for p in &lipid_particles[r.clone()] {
                sx += p.position[0];
                sy += p.position[1];
            }
            [sx / n, sy / n]
        })
        .collect();

    let mut disp = vec![[0.0_f64; 2]; n_mols];
    for _iter in 0..MAX_ITERS {
        disp.fill([0.0, 0.0]);
        let mut max_probe_overlap = 0.0_f64;

        // 1. Push any overlapping lipid rigidly in (x, y) away from probe beads.
        for (mi, r) in ranges.iter().enumerate() {
            for li in r.clone() {
                let lp = lipid_particles[li].position;
                let lr = bead_radii[li];
                for &(pp, pr) in probe {
                    let dz = lp[2] - pp[2];
                    let target = pr + lr + MARGIN;
                    if dz.abs() >= target {
                        continue;
                    }
                    let dx = pbc_delta(lp[0] - pp[0], bx);
                    let dy = pbc_delta(lp[1] - pp[1], by);
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    if dist < target {
                        let overlap = target - dist;
                        max_probe_overlap = max_probe_overlap.max(overlap);
                        let r_xy = dx.hypot(dy);
                        let (ux, uy) = if r_xy > 1e-6 {
                            (dx / r_xy, dy / r_xy)
                        } else {
                            let angle = (mi as f64 + 1.0) * 2.399_963_229_728_653;
                            (angle.cos(), angle.sin())
                        };
                        disp[mi][0] += 0.65 * overlap * ux;
                        disp[mi][1] += 0.65 * overlap * uy;
                    }
                }
            }
        }

        if max_probe_overlap < 0.01 {
            break;
        }

        // 2. Soft lateral separation between same-leaflet lipid neighbors so
        //    pushed lipids spread smoothly instead of piling up.
        for i in 0..n_mols {
            for j in (i + 1)..n_mols {
                if (mol_z_mean[i] - mol_z_mean[j]).abs() > 1.5 {
                    continue;
                }
                let dx = pbc_delta(mol_xy[i][0] - mol_xy[j][0], bx);
                let dy = pbc_delta(mol_xy[i][1] - mol_xy[j][1], by);
                let d_xy = dx.hypot(dy);
                const MIN_LIPID_SEP: f64 = 0.72;
                if d_xy < MIN_LIPID_SEP && d_xy > 1e-6 {
                    let push = 0.25 * (MIN_LIPID_SEP - d_xy);
                    let ux = dx / d_xy;
                    let uy = dy / d_xy;
                    disp[i][0] += push * ux;
                    disp[i][1] += push * uy;
                    disp[j][0] -= push * ux;
                    disp[j][1] -= push * uy;
                }
            }
        }

        // 3. Apply clamped lateral shifts rigidly to each lipid molecule.
        for (mi, r) in ranges.iter().enumerate() {
            let [mut sx, mut sy] = disp[mi];
            let mag = sx.hypot(sy);
            if mag < 1e-9 {
                continue;
            }
            if mag > MAX_STEP {
                let scale = MAX_STEP / mag;
                sx *= scale;
                sy *= scale;
            }
            mol_xy[mi][0] += sx;
            mol_xy[mi][1] += sy;
            for p in &mut lipid_particles[r.clone()] {
                p.position[0] += sx;
                p.position[1] += sy;
            }
        }
    }
}

/// Build the GENESIS `[SELECTION]` `group2` expression selecting all lipid
/// residues in the bilayer (e.g. `"rnam:POPC"` or `"rnam:POPC or rnam:POPS"`).
pub fn bilayer_selection(lipid_counts: &[(String, usize)]) -> String {
    let mut seen: Vec<&str> = Vec::new();
    for (name, _) in lipid_counts {
        if !seen.contains(&name.as_str()) {
            seen.push(name.as_str());
        }
    }
    seen.into_iter()
        .map(|name| format!("rnam:{name}"))
        .collect::<Vec<_>>()
        .join(" or ")
}

/// Auto-detected GENESIS umbrella output layout inside a campaign directory
/// (no `manifest.json` required).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CvSource {
    /// Independent per-window directories `win_1/<log>` .. `win_N/<log>`.
    Windows { count: usize },
    /// Multi-replica REUS logs `<prefix>_1.log` + `<prefix>_1.rem` .. `N`.
    Reus { count: usize },
}

/// Detect whether `dir` holds GENESIS REUS logs (`<reus_prefix>_1..N.{log,rem}`)
/// or independent-window logs (`win_1..N/<win_log>`).
pub fn detect_cv_source(dir: &Path, win_log: &str, reus_prefix: &str) -> Option<CvSource> {
    if dir.join(format!("{reus_prefix}_1.log")).is_file()
        && dir.join(format!("{reus_prefix}_1.rem")).is_file()
    {
        let mut count = 1usize;
        while dir
            .join(format!("{reus_prefix}_{}.log", count + 1))
            .is_file()
            && dir
                .join(format!("{reus_prefix}_{}.rem", count + 1))
                .is_file()
        {
            count += 1;
        }
        return Some(CvSource::Reus { count });
    }
    if dir.join("win_1").join(win_log).is_file() {
        let mut count = 1usize;
        while dir
            .join(format!("win_{}", count + 1))
            .join(win_log)
            .is_file()
        {
            count += 1;
        }
        return Some(CvSource::Windows { count });
    }
    None
}

/// Extract `(step, cv_value)` pairs from a GENESIS `atdyn` log file by locating
/// `STEP` and the first restraint CV column (`RESTR_CVS001`, `RESTR_CVS1`, …)
/// in the `INFO:` table header.
///
/// Rows with `step < discard_steps` are skipped.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] if the log contains no `INFO:` header
/// with `STEP` and `RESTR_CVS001` (or `RESTR_CVS1`).
pub fn extract_cv_from_log(
    log_content: &str,
    discard_steps: usize,
) -> Result<Vec<(usize, f64)>, BuildError> {
    let mut step_col: Option<usize> = None;
    let mut cv_col: Option<usize> = None;
    let mut rows = Vec::new();

    for line in log_content.lines() {
        let trimmed = line.trim_start();
        if !trimmed.starts_with("INFO:") {
            continue;
        }
        let toks: Vec<&str> = trimmed.split_whitespace().collect();
        if toks.len() < 3 {
            continue;
        }
        if toks[1] == "STEP" {
            step_col = Some(1);
            cv_col = toks.iter().position(|&t| {
                t.strip_prefix("RESTR_CVS")
                    .and_then(|num| num.parse::<usize>().ok())
                    == Some(1)
            });
            continue;
        }
        let (Some(sc), Some(cc)) = (step_col, cv_col) else {
            continue;
        };
        if toks.len() <= sc.max(cc) {
            continue;
        }
        let Ok(step) = toks[sc].parse::<usize>() else {
            continue;
        };
        if step < discard_steps {
            continue;
        }
        if let Ok(cv) = toks[cc].parse::<f64>()
            && cv.is_finite()
        {
            rows.push((step, cv));
        }
    }

    if cv_col.is_none() {
        return Err(BuildError::StructureParse(
            "GENESIS log contains no 'INFO: STEP ... RESTR_CVS001' header".to_string(),
        ));
    }
    Ok(rows)
}

/// Render extracted `(step, cv)` samples in GENESIS `mbar_analysis` `cvfile` format.
pub fn render_cv_dat(rows: &[(usize, f64)]) -> String {
    let mut s = String::with_capacity(rows.len() * 24);
    for &(step, cv) in rows {
        let _ = writeln!(s, "{step:>10} {cv:>14.6}");
    }
    s
}

/// Demultiplex GENESIS REUS `remfile` + `logfile` pairs into per-parameter-ID
/// CV series (`1..=n_replicas`) so `mbar_analysis` (`input_type = US`) can be
/// run directly on REUS output.
///
/// `rem_logs` is indexed by replica `0..n_replicas`; each entry is
/// `(log_content, rem_content)`.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] if any log is missing `RESTR_CVS001`.
pub fn demux_reus_cv(
    rem_logs: &[(&str, &str)],
    discard_steps: usize,
) -> Result<Vec<Vec<(usize, f64)>>, BuildError> {
    let n = rem_logs.len();
    let mut by_param: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

    for (rep_idx, &(log_txt, rem_txt)) in rem_logs.iter().enumerate() {
        let cvs = extract_cv_from_log(log_txt, discard_steps)?;
        // Parse `.rem` mapping: lines of `<step> <param_id>` (1-based).
        let mut param_map: Vec<(usize, usize)> = Vec::new();
        for line in rem_txt.lines() {
            let mut it = line.split_whitespace();
            if let (Some(st), Some(pid)) = (it.next(), it.next())
                && let (Ok(step), Ok(param_id)) = (st.parse::<usize>(), pid.parse::<usize>())
                && (1..=n).contains(&param_id)
            {
                param_map.push((step, param_id));
            }
        }
        // Track the active parameter ID for each logged step (initially rep_idx + 1).
        let mut map_pos = 0usize;
        let mut cur_param = rep_idx + 1;
        for (step, cv) in cvs {
            while map_pos < param_map.len() && param_map[map_pos].0 <= step {
                cur_param = param_map[map_pos].1;
                map_pos += 1;
            }
            by_param[cur_param - 1].push((step, cv));
        }
    }
    for series in &mut by_param {
        series.sort_unstable_by_key(|&(step, _)| step);
    }
    Ok(by_param)
}

/// Write-only summary metadata written to `manifest.json` so external Bayesian
/// optimization or HPC drivers can discover windows, paths, and restraint
/// parameters cleanly.
#[derive(Debug, Clone)]
pub struct UmbrellaManifest<'a> {
    pub probe: &'a str,
    pub probe_kind: ProbeKind,
    pub probe_atoms: usize,
    pub mode: &'a str,
    pub ensemble: &'a str,
    pub temperature_k: f64,
    pub k_kcal_mol_a2: f64,
    pub restraint_function: &'a str,
    pub restraint_direction: &'a str,
    pub group1: &'a str,
    pub group2: &'a str,
    pub box_nm: [f64; 3],
    pub lipid_counts: &'a [(String, usize)],
    pub z_windows_nm: &'a [f64],
}

impl UmbrellaManifest<'_> {
    /// Render pretty-printed JSON with zero external dependencies.
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(1024 + self.z_windows_nm.len() * 140);
        s.push_str("{\n");
        let _ = writeln!(s, "  \"probe\": \"{}\",", self.probe);
        let _ = writeln!(s, "  \"probe_kind\": \"{}\",", self.probe_kind.label());
        let _ = writeln!(s, "  \"probe_atoms\": {},", self.probe_atoms);
        let _ = writeln!(s, "  \"mode\": \"{}\",", self.mode);
        let _ = writeln!(s, "  \"ensemble\": \"{}\",", self.ensemble);
        let _ = writeln!(s, "  \"temperature_k\": {},", self.temperature_k);
        let _ = writeln!(s, "  \"k_kcal_mol_a2\": {},", self.k_kcal_mol_a2);
        let _ = writeln!(s, "  \"bias_convention\": \"U = k * (dz - z0)^2\",");
        let _ = writeln!(
            s,
            "  \"restraint_function\": \"{}\",",
            self.restraint_function
        );
        let _ = writeln!(
            s,
            "  \"restraint_direction\": \"{}\",",
            self.restraint_direction
        );
        let _ = writeln!(s, "  \"group1\": \"{}\",", self.group1);
        let _ = writeln!(s, "  \"group2\": \"{}\",", self.group2);
        let _ = writeln!(
            s,
            "  \"box_nm\": [{:.4}, {:.4}, {:.4}],",
            self.box_nm[0], self.box_nm[1], self.box_nm[2]
        );
        let mut aggregated: Vec<(&str, usize)> = Vec::new();
        for (name, count) in self.lipid_counts {
            if let Some((_, total)) = aggregated.iter_mut().find(|(n, _)| *n == name.as_str()) {
                *total += count;
            } else {
                aggregated.push((name.as_str(), *count));
            }
        }
        s.push_str("  \"lipids\": {");
        for (i, (name, count)) in aggregated.iter().enumerate() {
            if i > 0 {
                s.push_str(", ");
            }
            let _ = write!(s, "\"{name}\": {count}");
        }
        s.push_str("},\n");
        s.push_str("  \"topology\": {\n");
        s.push_str("    \"top\": \"system.top\",\n");
        s.push_str("    \"isolf_itp\": \"isolf.itp\",\n");
        s.push_str("    \"aicg2p_itp\": \"aicg2p.itp\",\n");
        s.push_str("    \"lipid_protein_itp\": \"lipid_protein.itp\",\n");
        let _ = writeln!(s, "    \"probe_itp\": \"{}.itp\"", self.probe);
        s.push_str("  },\n");
        s.push_str("  \"mbar_inp\": \"mbar.inp\",\n");
        s.push_str("  \"windows\": [\n");
        let is_reus = self.mode == "reus";
        for (idx, &z_nm) in self.z_windows_nm.iter().enumerate() {
            let id = idx + 1;
            let z_a = z_nm * 10.0;
            let comma = if id < self.z_windows_nm.len() {
                ","
            } else {
                ""
            };
            if is_reus {
                let _ = writeln!(
                    s,
                    "    {{\"id\": {id}, \"z_nm\": {z_nm:.4}, \"z_angstrom\": {z_a:.4}, \
                     \"gro\": \"win_{id}.gro\", \"min_inp\": \"min_{id}.inp\", \
                     \"log\": \"prod_reus_{id}.log\", \"rem\": \"prod_reus_{id}.rem\", \
                     \"cv_dat\": \"cv/{id}.dat\"}}{comma}"
                );
            } else {
                let _ = writeln!(
                    s,
                    "    {{\"id\": {id}, \"z_nm\": {z_nm:.4}, \"z_angstrom\": {z_a:.4}, \
                     \"dir\": \"win_{id}\", \"gro\": \"win_{id}/system.gro\", \
                     \"min_inp\": \"win_{id}/min.inp\", \"equil_inp\": \"win_{id}/equil.inp\", \
                     \"prod_inp\": \"win_{id}/prod.inp\", \"log\": \"win_{id}/prod.log\", \
                     \"cv_dat\": \"cv/{id}.dat\"}}{comma}"
                );
            }
        }
        s.push_str("  ]\n");
        s.push_str("}\n");
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composition::Composition;
    use crate::membrane::{BuildOptions, Sizing};
    use crate::protein::aicg2p::ProteinOptions;
    use crate::protein::structure::parse_pdb;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn sample_helix_protein() -> Protein {
        let mut pdb = String::new();
        let mut serial = 1;
        let seq = [
            "ALA", "LEU", "VAL", "ILE", "LEU", "PHE", "TRP", "LEU", "LYS", "ALA",
        ];
        for (i, &res) in seq.iter().enumerate() {
            let idx = i + 1;
            let theta = i as f64 * (100.0_f64).to_radians();
            let r = 2.3;
            let z = (i as f64 - 4.5) * 1.5;
            let ca = [r * theta.cos(), r * theta.sin(), z];
            let n = [r * (theta - 0.3).cos(), r * (theta - 0.3).sin(), z - 0.6];
            let c = [r * (theta + 0.3).cos(), r * (theta + 0.3).sin(), z + 0.6];
            for (aname, p) in [("N", n), ("CA", ca), ("C", c)] {
                let _ = writeln!(
                    pdb,
                    "ATOM  {serial:>5}  {aname:<3} {res:>3} A{idx:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00           C",
                    p[0], p[1], p[2]
                );
                serial += 1;
            }
        }
        pdb.push_str("END\n");
        let structure = parse_pdb(&pdb).unwrap();
        Protein::from_structure("helix", &structure, &ProteinOptions::default()).unwrap()
    }

    #[test]
    fn parses_z_window_ranges_counts_and_lists() {
        let w = parse_z_windows("0.0..4.0:0.1").unwrap();
        assert_eq!(w.len(), 41);
        assert_eq!(w[0], 0.0);
        assert_eq!(w[10], 1.0);
        assert_eq!(w[40], 4.0);

        let w2 = parse_z_windows("-2.0..2.0/5").unwrap();
        assert_eq!(w2, vec![-2.0, -1.0, 0.0, 1.0, 2.0]);

        let w3 = parse_z_windows("0, 0.5, 1.25, 3.5").unwrap();
        assert_eq!(w3, vec![0.0, 0.5, 1.25, 3.5]);

        assert!(parse_z_windows("2.0..1.0:0.1").is_err());
        assert!(parse_z_windows("0.0..1.0:0").is_err());
        assert!(parse_z_windows("0.0..1.0/0").is_err());
    }

    #[test]
    fn centers_pbc_split_bilayer_intact_at_origin() {
        let ff = ForceField::isolf();
        let comp = Composition::parse("POPC").unwrap();
        let mut rng = StdRng::seed_from_u64(4);
        let mut mem = Membrane::build(
            &ff,
            "popc",
            &comp,
            &comp,
            Sizing::Count {
                upper: 16,
                lower: 16,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap();
        mem.box_size[2] = 20.0;
        // Shift bilayer so it wraps across the z = 0 / z = 20 periodic boundary.
        for p in &mut mem.particles {
            p.position[2] = (p.position[2] + 20.0) % 20.0;
        }
        center_bilayer_at_origin(&mut mem, 20.0);
        let z_min = mem
            .particles
            .iter()
            .map(|p| p.position[2])
            .fold(f64::INFINITY, f64::min);
        let z_max = mem
            .particles
            .iter()
            .map(|p| p.position[2])
            .fold(f64::NEG_INFINITY, f64::max);
        // Bilayer slab (~4 nm thick) must stay contiguous around z = 0, not split across ±10 nm.
        assert!(
            z_max - z_min < 6.0,
            "expected contiguous slab, got [{z_min}, {z_max}]"
        );
        let (msum, zsum) = mem.particles.iter().fold((0.0, 0.0), |(m, z), p| {
            (m + p.mass, z + p.mass * p.position[2])
        });
        assert!((zsum / msum).abs() < 1e-10);
    }

    #[test]
    fn window_builder_works_for_single_aa_and_multi_residue_peptide_without_carving() {
        let ff = ForceField::isolf();
        let comp = Composition::parse("POPC").unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        let mut mem = Membrane::build(
            &ff,
            "popc",
            &comp,
            &comp,
            Sizing::Count {
                upper: 32,
                lower: 32,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap();
        center_bilayer_at_origin(&mut mem, 20.0);

        // 1. Single amino-acid bead.
        let trp = Protein::single_bead("TRP", None).unwrap();
        let parts = build_window_particles(&mem, &trp, 1.5, &ff);
        assert_eq!(parts.len(), 1 + mem.particles.len());
        assert_eq!(parts[0].residue, "TRP");
        assert_eq!(probe_selection(trp.beads.len()), "ano:1");
        assert!((parts[0].position[2] - 1.5).abs() < 1e-12);

        // 2. Multi-residue helix at z = 0.0 (bilayer midplane): zero lipids carved,
        //    exact probe COM at z = 0, bilayer COM still at z = 0, lipids declashed.
        let helix = sample_helix_protein();
        let n_p = helix.beads.len();
        assert_eq!(probe_selection(n_p), format!("ano:1-{n_p}"));
        let sys = build_window_system("helix win".to_string(), &mem, &helix, 0.0, &ff);
        assert_eq!(sys.particles.len(), n_p + mem.particles.len());
        assert_eq!(sys.total_lipids(), mem.total_lipids());
        let (pm, pz) = sys.particles[..n_p].iter().fold((0.0, 0.0), |(m, z), p| {
            (m + p.mass, z + p.mass * p.position[2])
        });
        assert!((pz / pm).abs() < 1e-10);
        let (lm, lz) = sys.particles[n_p..].iter().fold((0.0, 0.0), |(m, z), p| {
            (m + p.mass, z + p.mass * p.position[2])
        });
        assert!((lz / lm).abs() < 1e-10);
        let top = sys.to_top_with_lipid_protein(true);
        assert!(top.contains("#include \"./lipid_protein.itp\""));
        assert!(top.contains("#include \"./helix.itp\""));
    }

    #[test]
    fn extracts_restr_cvs_demuxes_reus_and_renders_manifest() {
        let log1 = "\
INFO:       STEP            TIME       TOTAL_ENE   RESTRAINT_TOTAL   RESTRAINT001     RESTR_CVS001
------------ --------------- --------------- --------------- --------------- ---------------
INFO:          0          0.0000        -100.000             0.010             0.010           10.1414
INFO:        100          1.0000        -101.000             0.020             0.020            0.1000
INFO:        200          2.0000        -102.000             0.005             0.005            5.1000
";
        let log2 = "\
INFO: STEP RESTR_CVS001
INFO: 100  4.90
INFO: 200 -0.10
";
        let rows = extract_cv_from_log(log1, 100).unwrap();
        assert_eq!(rows, vec![(100, 0.10), (200, 5.10)]);
        assert!(render_cv_dat(&rows).contains("100       0.100000"));

        let rem1 = "100 1\n200 2\n";
        let rem2 = "100 2\n200 1\n";
        let demux = demux_reus_cv(&[(log1, rem1), (log2, rem2)], 100).unwrap();
        assert_eq!(demux[0], vec![(100, 0.10), (200, -0.10)]);
        assert_eq!(demux[1], vec![(100, 4.90), (200, 5.10)]);

        let counts = vec![
            ("POPC".to_string(), 32),
            ("POPS".to_string(), 16),
            ("POPC".to_string(), 32),
        ];
        assert_eq!(bilayer_selection(&counts), "rnam:POPC or rnam:POPS");
        let manifest = UmbrellaManifest {
            probe: "TRP",
            probe_kind: ProbeKind::AminoAcid,
            probe_atoms: 1,
            mode: "reus",
            ensemble: "nvt",
            temperature_k: 323.15,
            k_kcal_mol_a2: 0.5,
            restraint_function: "DISTMASS",
            restraint_direction: "Z2",
            group1: "ano:1",
            group2: "rnam:POPC or rnam:POPS",
            box_nm: [7.0, 7.0, 20.0],
            lipid_counts: &counts,
            z_windows_nm: &[0.0, 0.5, 1.0],
        };
        let json = manifest.to_json();
        assert!(json.contains("\"probe_kind\": \"amino_acid\""));
        assert!(json.contains("\"group1\": \"ano:1\""));
        assert!(json.contains("\"lipids\": {\"POPC\": 64, \"POPS\": 16}"));
    }
}
