//! AICG2+ and Clementi Gō protein coarse-graining and GENESIS `.itp` rendering.
//!
//! Faithfully ports `genesis_cg_tool`'s protein builder (`coarse_graining.jl`,
//! `coarse_graining_subroutines.jl`, `biomath.jl`, `parser_top.jl`):
//!
//! - Heavy-atom classification (17 AICG2+ pairwise contact classes with
//!   salt-bridge cap and energy clamping to `[-5.0, -0.5]` kcal/mol).
//! - Exact 3D spatial-grid + residue bounding-sphere native contact search
//!   (`d < 6.5 Å`), parallelized with Rayon and deterministically ordered.
//! - General (`aicg-scale = 1`, default) and average (`aicg-scale = 0`) energy
//!   scaling, flexible-local angle/dihedral terms, safe dihedral modes
//!   (`0`/`1`/`3`), IDR local/nonlocal exclusions, and RESPAC charge overrides.

use std::collections::HashSet;
use std::fmt::Write;

use rayon::prelude::*;

use crate::error::BuildError;
use crate::protein::params::{
    AICG_13_AVE, AICG_13_GEN, AICG_13_SIGMA, AICG_14_AVE, AICG_14_GEN, AICG_14_SIGMA,
    AICG_ATOMIC_CUTOFF, AICG_BOND_K, AICG_CONTACT_AVE, AICG_CONTACT_GEN, AICG_ENE_LOWER_LIM,
    AICG_ENE_UPPER_LIM, AICG_GO_ATOMIC_CUTOFF, AICG_HYDROGEN_BOND_CUTOFF, AICG_PAIRWISE_ENERGY,
    AICG_SALT_BRIDGE_CUTOFF, CAL_TO_KJ, CCGO_ANGL_K, CCGO_BOND_K, CCGO_DIHE_K_1, CCGO_DIHE_K_3,
    CCGO_NATIVE_EPSILON, ITYPE_BB_CX, ITYPE_BB_DA, ITYPE_BB_HB, ITYPE_BB_XX, ITYPE_LR_CT,
    ITYPE_OFFST, ITYPE_SB_CX, ITYPE_SB_DA, ITYPE_SB_HB, ITYPE_SB_QX, ITYPE_SB_XX, ITYPE_SS_CX,
    ITYPE_SS_DA, ITYPE_SS_HB, ITYPE_SS_QX, ITYPE_SS_SB, ITYPE_SS_XX, aa_charge, aa_mass, aa_rmin,
    fmt_sci,
};
use crate::protein::structure::{AllAtom, AllAtomStructure};
use crate::protein::{
    Protein, ProteinAngle, ProteinBead, ProteinBond, ProteinChain, ProteinDihedral, ProteinPair,
};

/// Protein coarse-grained force-field model.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ProteinModel {
    /// AICG2+ structure-based model with flexible local potentials (default).
    #[default]
    Aicg2Plus,
    /// Clementi Gō-like Cα model with harmonic angles and periodic dihedrals.
    ClementiGo,
}

/// Energy scaling scheme for AICG2+ local and contact energies.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum AicgScale {
    /// Universal scaling factors (`1.11`, `0.87`, `0.32`; `aicg-scale = 1`, default).
    #[default]
    General,
    /// Normalize mean energies to reference averages (`1.72`, `1.23`, `0.55`; `aicg-scale = 0`).
    Average,
}

/// Safe dihedral potential damping mode (matches `genesis_cg_tool --use-safe-dihedral`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SafeDihedral {
    /// Unmodified dihedral potentials (`0`: Gaussian `21`, tabulated `22`, periodic `1`).
    None,
    /// `cos²(kθ)`-damped safe dihedrals (`1`, default: Gaussian `41`, tabulated `52`, periodic `32`).
    #[default]
    Cos2,
    /// `sin³(kθ)`-damped safe dihedrals (`3`: Gaussian `43`, tabulated `22`, periodic `33`).
    Sin3,
}

/// Options controlling protein coarse-graining.
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinOptions {
    /// Protein force field (`AICG2+` or `Clementi`).
    pub model: ProteinModel,
    /// AICG2+ energy scaling scheme (`General` or `Average`).
    pub scale: AicgScale,
    /// Dihedral singularity-damping mode.
    pub safe_dihedral: SafeDihedral,
    /// 1-based bead indices excluded from AICG2+ local Gaussian 1–3 and 1–4 terms.
    pub idr_local: Vec<usize>,
    /// 1-based bead indices excluded from native Gō contacts (`[ pairs ]`).
    pub idr_nonlocal: Vec<usize>,
    /// Optional `(1-based bead index, charge)` overrides (e.g. RESPAC charges).
    pub charges: Vec<(usize, f64)>,
    /// Scaling multiplier for Clementi Gō native contact strength (default `1.0`).
    pub ccgo_contact_scale: f64,
}

impl Default for ProteinOptions {
    fn default() -> Self {
        Self {
            model: ProteinModel::Aicg2Plus,
            scale: AicgScale::General,
            safe_dihedral: SafeDihedral::Cos2,
            idr_local: Vec::new(),
            idr_nonlocal: Vec::new(),
            charges: Vec::new(),
            ccgo_contact_scale: 1.0,
        }
    }
}

/// Parse a residue/bead index selection such as `"1-25,40-55"`, `"1 to 25, 40"`,
/// or `"10"` into a sorted, deduplicated list of 1-based indices.
///
/// # Errors
///
/// Returns [`BuildError::InvalidSelection`] on non-positive numbers or reversed ranges.
pub fn parse_selection(spec: &str) -> Result<Vec<usize>, BuildError> {
    let mut out = Vec::new();
    for raw_part in spec.split(',') {
        let part = raw_part.trim();
        if part.is_empty() {
            continue;
        }
        let split_pair = if let Some((a, b)) = part.split_once(" to ") {
            Some((a, b))
        } else if let Some((a, b)) = part.split_once(" TO ") {
            Some((a, b))
        } else {
            part.split_once('-')
        };
        if let Some((a, b)) = split_pair {
            let start = parse_pos_index(a.trim(), part)?;
            let end = parse_pos_index(b.trim(), part)?;
            if end < start {
                return Err(BuildError::InvalidSelection(format!(
                    "reversed range '{part}' ({start} > {end})"
                )));
            }
            out.extend(start..=end);
        } else {
            out.push(parse_pos_index(part, part)?);
        }
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

fn parse_pos_index(s: &str, context: &str) -> Result<usize, BuildError> {
    s.parse::<usize>().ok().filter(|&v| v > 0).ok_or_else(|| {
        BuildError::InvalidSelection(format!("'{context}' is not a valid 1-based index/range"))
    })
}

/// Parse a whitespace-delimited charge override table (`<1-based idx> <charge>` per line).
/// Blank lines and comments starting with `#` or `;` are ignored.
///
/// # Errors
///
/// Returns [`BuildError::InvalidSelection`] on malformed lines.
pub fn parse_charge_table(content: &str) -> Result<Vec<(usize, f64)>, BuildError> {
    let mut out = Vec::new();
    for (line_no, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with(';') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let idx_str = parts.next().unwrap_or("");
        let q_str = parts.next().ok_or_else(|| {
            BuildError::InvalidSelection(format!(
                "charge file line {}: expected '<index> <charge>'",
                line_no + 1
            ))
        })?;
        let idx = idx_str
            .parse::<usize>()
            .ok()
            .filter(|&v| v > 0)
            .ok_or_else(|| {
                BuildError::InvalidSelection(format!(
                    "charge file line {}: invalid 1-based bead index '{idx_str}'",
                    line_no + 1
                ))
            })?;
        let q = q_str
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| {
                BuildError::InvalidSelection(format!(
                    "charge file line {}: invalid charge '{q_str}'",
                    line_no + 1
                ))
            })?;
        out.push((idx, q));
    }
    Ok(out)
}

// ==================== Geometry helpers (Å) ====================

#[inline]
fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

#[inline]
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    dist_sq(a, b).sqrt()
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Angle `p1–p2–p3` in degrees (matches `genesis_cg_tool/src/lib/biomath.jl`).
fn compute_angle_deg(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let v1 = sub(p1, p2);
    let v2 = sub(p3, p2);
    let n1 = dot(v1, v1).sqrt();
    let n2 = dot(v2, v2).sqrt();
    if n1 == 0.0 || n2 == 0.0 {
        return 0.0;
    }
    (dot(v1, v2) / n1 / n2).clamp(-1.0, 1.0).acos().to_degrees()
}

/// Dihedral `p1–p2–p3–p4` in degrees (matches `genesis_cg_tool/src/lib/biomath.jl`).
fn compute_dihedral_deg(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let v12 = sub(p2, p1);
    let v23 = sub(p3, p2);
    let v34 = sub(p4, p3);
    let c123 = cross(v12, v23);
    let c234 = cross(v23, v34);
    let n1 = dot(c123, c123).sqrt();
    let n2 = dot(c234, c234).sqrt();
    if n1 == 0.0 || n2 == 0.0 {
        return 0.0;
    }
    let dih = (dot(c123, c234) / n1 / n2)
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees();
    let judge = dot(cross(c123, c234), v23);
    if judge < 0.0 { -dih } else { dih }
}

// ==================== AICG2+ heavy-atom classification ====================

#[derive(Clone, Copy)]
struct HeavyAtom {
    pos: [f64; 3],
    is_backbone: bool,
    is_donor: bool,
    is_acceptor: bool,
    is_cation: bool,
    is_anion: bool,
    is_carbon: bool,
}

impl HeavyAtom {
    fn from_all_atom(atom: &AllAtom, res_name: &str) -> Option<Self> {
        let first = atom.name.as_bytes().first().copied()?;
        if first == b'H' {
            return None;
        }
        let name = atom.name.as_str();
        let is_backbone = matches!(name, "N" | "C" | "O" | "OXT" | "CA");
        let is_donor = first == b'N'
            || (first == b'S' && res_name == "CYS")
            || (first == b'O'
                && ((res_name == "SER" && name == "OG")
                    || (res_name == "THR" && name == "OG1")
                    || (res_name == "TYR" && name == "OH")));
        let is_acceptor = first == b'O' || first == b'S';
        let is_cation = first == b'N'
            && ((res_name == "ARG" && matches!(name, "NH1" | "NH2"))
                || (res_name == "LYS" && name == "NZ"));
        let is_anion = first == b'O'
            && ((res_name == "GLU" && matches!(name, "OE1" | "OE2"))
                || (res_name == "ASP" && matches!(name, "OD1" | "OD2")));
        let is_carbon = first == b'C';
        Some(Self {
            pos: atom.pos,
            is_backbone,
            is_donor,
            is_acceptor,
            is_cation,
            is_anion,
            is_carbon,
        })
    }
}

/// Precomputed heavy-atom list and bounding sphere (Å) for one coarse-grained residue.
struct ResidueHeavy {
    atoms: Vec<HeavyAtom>,
    center: [f64; 3],
    radius: f64,
    ca_pos: [f64; 3],
    chain_idx: usize,
}

/// True if any heavy-atom pair between `a` and `b` is closer than `6.5 Å`.
fn is_native_contact(a: &ResidueHeavy, b: &ResidueHeavy) -> bool {
    let bound = a.radius + b.radius + AICG_GO_ATOMIC_CUTOFF;
    if dist_sq(a.center, b.center) >= bound * bound {
        return false;
    }
    let cut2 = AICG_GO_ATOMIC_CUTOFF * AICG_GO_ATOMIC_CUTOFF;
    for ha in &a.atoms {
        for hb in &b.atoms {
            if dist_sq(ha.pos, hb.pos) < cut2 {
                return true;
            }
        }
    }
    false
}

/// Compute the raw clamped AICG2+ pairwise energy (kcal/mol) between two
/// residues (`[-5.0, -0.5]`), matching `count_aicg_atomic_contact` verbatim.
fn aicg_raw_energy(a: &[HeavyAtom], b: &[HeavyAtom]) -> f64 {
    let mut counts = [0i32; 17];
    counts[ITYPE_OFFST] = 1;
    let mut short_range = 0i32;

    let go_cut2 = AICG_GO_ATOMIC_CUTOFF * AICG_GO_ATOMIC_CUTOFF;
    for ha in a {
        for hb in b {
            let d2 = dist_sq(ha.pos, hb.pos);
            if d2 >= go_cut2 {
                continue;
            }
            counts[ITYPE_LR_CT] += 1;
            let d = d2.sqrt();
            if d < AICG_ATOMIC_CUTOFF {
                short_range += 1;
                let is_hb = (ha.is_acceptor && hb.is_donor) || (hb.is_acceptor && ha.is_donor);
                let is_sb = (ha.is_cation && hb.is_anion) || (hb.is_cation && ha.is_anion);
                let is_nonsb_charge = ha.is_cation || ha.is_anion || hb.is_cation || hb.is_anion;

                if ha.is_backbone && hb.is_backbone {
                    if is_hb {
                        if d < AICG_HYDROGEN_BOND_CUTOFF {
                            counts[ITYPE_BB_HB] += 1;
                        } else {
                            counts[ITYPE_BB_DA] += 1;
                        }
                    } else if ha.is_carbon || hb.is_carbon {
                        counts[ITYPE_BB_CX] += 1;
                    } else {
                        counts[ITYPE_BB_XX] += 1;
                    }
                } else if !ha.is_backbone && !hb.is_backbone {
                    if is_hb {
                        if is_sb {
                            if d < AICG_SALT_BRIDGE_CUTOFF {
                                counts[ITYPE_SS_SB] += 1;
                            } else {
                                counts[ITYPE_SS_QX] += 1;
                            }
                        } else if d < AICG_HYDROGEN_BOND_CUTOFF {
                            counts[ITYPE_SS_HB] += 1;
                        } else if is_nonsb_charge {
                            counts[ITYPE_SS_QX] += 1;
                        } else {
                            counts[ITYPE_SS_DA] += 1;
                        }
                    } else if is_nonsb_charge {
                        counts[ITYPE_SS_QX] += 1;
                    } else if ha.is_carbon || hb.is_carbon {
                        counts[ITYPE_SS_CX] += 1;
                    } else {
                        counts[ITYPE_SS_XX] += 1;
                    }
                } else if is_hb {
                    if d < AICG_HYDROGEN_BOND_CUTOFF {
                        counts[ITYPE_SB_HB] += 1;
                    } else if is_nonsb_charge {
                        counts[ITYPE_SB_QX] += 1;
                    } else {
                        counts[ITYPE_SB_DA] += 1;
                    }
                } else if is_nonsb_charge {
                    counts[ITYPE_SB_QX] += 1;
                } else if ha.is_carbon || hb.is_carbon {
                    counts[ITYPE_SB_CX] += 1;
                } else {
                    counts[ITYPE_SB_XX] += 1;
                }
            }
        }
    }

    counts[ITYPE_LR_CT] -= short_range;
    if counts[ITYPE_SS_SB] >= 2 {
        counts[ITYPE_SS_QX] += counts[ITYPE_SS_SB] - 1;
        counts[ITYPE_SS_SB] = 1;
    }

    let mut energy = 0.0_f64;
    for (w, &c) in AICG_PAIRWISE_ENERGY.iter().zip(&counts) {
        energy += w * c as f64;
    }
    energy.clamp(AICG_ENE_LOWER_LIM, AICG_ENE_UPPER_LIM)
}

// ==================== Main builder ====================

/// Build a coarse-grained [`Protein`] (AICG2+ or Clementi Gō) from an
/// all-atom structure.
///
/// # Errors
///
/// Returns [`BuildError::EmptyProtein`] if there are no residues,
/// [`BuildError::StructureParse`] if a residue is missing a `CA` atom, or
/// [`BuildError::InvalidSelection`] if a charge override index is out of range.
pub fn build_protein(
    name: &str,
    structure: &AllAtomStructure,
    options: &ProteinOptions,
) -> Result<Protein, BuildError> {
    if structure.chains.is_empty() {
        return Err(BuildError::EmptyProtein);
    }

    let mut beads = Vec::new();
    let mut chains = Vec::new();
    let mut heavy_res = Vec::new();

    for (chain_idx, aa_chain) in structure.chains.iter().enumerate() {
        let start = beads.len();
        for &res_idx in &aa_chain.residues {
            let res = &structure.residues[res_idx];
            let mut ca_pos = None;
            let mut heavy = Vec::with_capacity(res.atoms.len());
            for &atom_idx in &res.atoms {
                let atom = &structure.atoms[atom_idx];
                if atom.name == "CA" && ca_pos.is_none() {
                    ca_pos = Some(atom.pos);
                }
                if let Some(ha) = HeavyAtom::from_all_atom(atom, &res.name) {
                    heavy.push(ha);
                }
            }
            let ca = ca_pos.ok_or_else(|| {
                BuildError::StructureParse(format!(
                    "chain '{}' residue {} {} has no CA atom",
                    aa_chain.chain_id, res.seq_id, res.name
                ))
            })?;

            // Bounding sphere of heavy atoms for exact contact acceleration.
            let mut center = ca;
            if !heavy.is_empty() {
                let inv = 1.0 / heavy.len() as f64;
                let mut sum = [0.0; 3];
                for ha in &heavy {
                    sum[0] += ha.pos[0];
                    sum[1] += ha.pos[1];
                    sum[2] += ha.pos[2];
                }
                center = [sum[0] * inv, sum[1] * inv, sum[2] * inv];
            }
            let radius = heavy
                .iter()
                .map(|ha| dist_sq(center, ha.pos))
                .fold(0.0_f64, f64::max)
                .sqrt();

            let bead_idx = beads.len() + 1;
            let rmin = aa_rmin(&res.name);
            beads.push(ProteinBead {
                residue: res.name.clone(),
                atom: "CA".to_string(),
                atom_type: res.name.clone(),
                res_id: bead_idx,
                atom_id: bead_idx,
                chain_idx,
                chain_id: aa_chain.chain_id,
                seg_name: aa_chain.seg_name.clone(),
                // Store coordinates in nm (PDB Å × 0.1) to match all of `isolf`.
                position: [ca[0] * 0.1, ca[1] * 0.1, ca[2] * 0.1],
                charge: aa_charge(&res.name),
                mass: aa_mass(&res.name),
                radius: 0.5 * rmin,
            });
            heavy_res.push(ResidueHeavy {
                atoms: heavy,
                center,
                radius,
                ca_pos: ca,
                chain_idx,
            });
        }
        chains.push(ProteinChain {
            chain_id: aa_chain.chain_id,
            seg_name: aa_chain.seg_name.clone(),
            range: start..beads.len(),
        });
    }

    if beads.is_empty() {
        return Err(BuildError::EmptyProtein);
    }

    // Apply optional RESPAC / user charge overrides (1-based indices).
    for &(idx, q) in &options.charges {
        if idx == 0 || idx > beads.len() {
            return Err(BuildError::InvalidSelection(format!(
                "charge override index {idx} is out of range 1..={}",
                beads.len()
            )));
        }
        beads[idx - 1].charge = q;
    }

    let idr_local: HashSet<usize> = options.idr_local.iter().copied().collect();
    let idr_nonlocal: HashSet<usize> = options.idr_nonlocal.iter().copied().collect();

    // ---- Bonds ----
    let mut bonds = Vec::new();
    let mut warnings = Vec::new();
    let bond_k = match options.model {
        ProteinModel::Aicg2Plus => AICG_BOND_K,
        ProteinModel::ClementiGo => CCGO_BOND_K,
    };
    for (ci, ch) in chains.iter().enumerate() {
        if ch.range.len() < 2 {
            continue;
        }
        for i in ch.range.start..ch.range.end - 1 {
            let r0 = dist(heavy_res[i].ca_pos, heavy_res[i + 1].ca_pos);
            if !(3.6..=4.0).contains(&r0) {
                warnings.push(format!(
                    "abnormal CA-CA bond in chain {} between beads {} and {}: {:.3} nm",
                    ci + 1,
                    i + 1,
                    i + 2,
                    r0 * 0.1
                ));
            }
            bonds.push(ProteinBond {
                i: i + 1,
                j: i + 2,
                func: 1,
                eq_nm: r0 * 0.1,
                coef_kj: bond_k * CAL_TO_KJ * 200.0,
            });
        }
    }

    // ---- Angles (1–3) ----
    let mut raw_13: Vec<(usize, usize, usize, f64, f64, f64)> = Vec::new();
    for ch in &chains {
        if ch.range.len() < 3 {
            continue;
        }
        for i in ch.range.start..ch.range.end - 2 {
            let d13 = dist(heavy_res[i].ca_pos, heavy_res[i + 2].ca_pos);
            let ang = compute_angle_deg(
                heavy_res[i].ca_pos,
                heavy_res[i + 1].ca_pos,
                heavy_res[i + 2].ca_pos,
            );
            let e13 = if options.model == ProteinModel::Aicg2Plus {
                aicg_raw_energy(&heavy_res[i].atoms, &heavy_res[i + 2].atoms)
            } else {
                0.0
            };
            raw_13.push((i + 1, i + 2, i + 3, d13, ang, e13));
        }
    }

    // ---- Dihedrals (1–4) ----
    let mut raw_14: Vec<(usize, usize, usize, usize, f64, f64)> = Vec::new();
    for ch in &chains {
        if ch.range.len() < 4 {
            continue;
        }
        for i in ch.range.start..ch.range.end - 3 {
            let dih = compute_dihedral_deg(
                heavy_res[i].ca_pos,
                heavy_res[i + 1].ca_pos,
                heavy_res[i + 2].ca_pos,
                heavy_res[i + 3].ca_pos,
            );
            let e14 = if options.model == ProteinModel::Aicg2Plus {
                aicg_raw_energy(&heavy_res[i].atoms, &heavy_res[i + 3].atoms)
            } else {
                0.0
            };
            raw_14.push((i + 1, i + 2, i + 3, i + 4, dih, e14));
        }
    }

    // Normalize AICG2+ local energies (1–3 and 1–4).
    if options.model == ProteinModel::Aicg2Plus {
        match options.scale {
            AicgScale::General => {
                for entry in &mut raw_13 {
                    entry.5 *= -AICG_13_GEN;
                }
                for entry in &mut raw_14 {
                    entry.5 *= -AICG_14_GEN;
                }
            }
            AicgScale::Average => {
                if !raw_13.is_empty() {
                    let mean13 = raw_13.iter().map(|e| e.5).sum::<f64>() / raw_13.len() as f64;
                    if mean13 != 0.0 {
                        let factor = AICG_13_AVE / mean13;
                        for entry in &mut raw_13 {
                            entry.5 *= factor;
                        }
                    }
                }
                if !raw_14.is_empty() {
                    let mean14 = raw_14.iter().map(|e| e.5).sum::<f64>() / raw_14.len() as f64;
                    if mean14 != 0.0 {
                        let factor = AICG_14_AVE / mean14;
                        for entry in &mut raw_14 {
                            entry.5 *= factor;
                        }
                    }
                }
            }
        }
    }

    let mut angles = Vec::new();
    let mut dihedrals = Vec::new();
    match options.model {
        ProteinModel::Aicg2Plus => {
            // 1. Gaussian 1–3 angles (f = 21), skipping IDR-local residues.
            for &(i, j, k, d13, _, e13) in &raw_13 {
                if idr_local.contains(&i) || idr_local.contains(&j) || idr_local.contains(&k) {
                    continue;
                }
                angles.push(ProteinAngle {
                    i,
                    j,
                    k,
                    func: 21,
                    eq: d13 * 0.1,
                    coef_kj: e13 * CAL_TO_KJ,
                    width: AICG_13_SIGMA * 0.1,
                });
            }
            // 2. Flexible-local angles (f = 22) for every consecutive triple.
            for &(i, j, k, _, _, _) in &raw_13 {
                angles.push(ProteinAngle {
                    i,
                    j,
                    k,
                    func: 22,
                    eq: 0.0,
                    coef_kj: 0.0,
                    width: 0.0,
                });
            }

            let gaus_func = match options.safe_dihedral {
                SafeDihedral::None => 21,
                SafeDihedral::Cos2 => 41,
                SafeDihedral::Sin3 => 43,
            };
            let tabu_func = match options.safe_dihedral {
                SafeDihedral::None | SafeDihedral::Sin3 => 22,
                SafeDihedral::Cos2 => 52,
            };

            // 1. Gaussian 1–4 dihedrals, skipping IDR-local residues.
            for &(i, j, k, l, dih, e14) in &raw_14 {
                if idr_local.contains(&i)
                    || idr_local.contains(&j)
                    || idr_local.contains(&k)
                    || idr_local.contains(&l)
                {
                    continue;
                }
                dihedrals.push(ProteinDihedral {
                    i,
                    j,
                    k,
                    l,
                    func: gaus_func,
                    eq_deg: dih,
                    coef_kj: e14 * CAL_TO_KJ,
                    width: AICG_14_SIGMA,
                    multiplicity: 0,
                });
            }
            // 2. Flexible-local dihedrals for every consecutive quadruple.
            for &(i, j, k, l, _, _) in &raw_14 {
                dihedrals.push(ProteinDihedral {
                    i,
                    j,
                    k,
                    l,
                    func: tabu_func,
                    eq_deg: 0.0,
                    coef_kj: 0.0,
                    width: 0.0,
                    multiplicity: 0,
                });
            }
        }
        ProteinModel::ClementiGo => {
            for &(i, j, k, _, ang, _) in &raw_13 {
                angles.push(ProteinAngle {
                    i,
                    j,
                    k,
                    func: 1,
                    eq: ang,
                    coef_kj: CCGO_ANGL_K * 2.0 * CAL_TO_KJ,
                    width: 0.0,
                });
            }
            let peri_func = match options.safe_dihedral {
                SafeDihedral::None => 1,
                SafeDihedral::Cos2 => 32,
                SafeDihedral::Sin3 => 33,
            };
            for &(i, j, k, l, dih, _) in &raw_14 {
                dihedrals.push(ProteinDihedral {
                    i,
                    j,
                    k,
                    l,
                    func: peri_func,
                    eq_deg: dih - 180.0,
                    coef_kj: CCGO_DIHE_K_1 * CAL_TO_KJ,
                    width: 0.0,
                    multiplicity: 1,
                });
            }
            for &(i, j, k, l, dih, _) in &raw_14 {
                dihedrals.push(ProteinDihedral {
                    i,
                    j,
                    k,
                    l,
                    func: peri_func,
                    eq_deg: 3.0 * dih - 180.0,
                    coef_kj: CCGO_DIHE_K_3 * CAL_TO_KJ,
                    width: 0.0,
                    multiplicity: 3,
                });
            }
        }
    }

    // ---- Native Gō contacts (intra-chain j >= i+4 and inter-chain i < j) ----
    let mut raw_contacts = find_native_contacts(&heavy_res, options.model);

    if options.model == ProteinModel::Aicg2Plus {
        match options.scale {
            AicgScale::General => {
                for c in &mut raw_contacts {
                    c.3 *= -AICG_CONTACT_GEN;
                }
            }
            AicgScale::Average => {
                if !raw_contacts.is_empty() {
                    let mean_c =
                        raw_contacts.iter().map(|c| c.3).sum::<f64>() / raw_contacts.len() as f64;
                    if mean_c != 0.0 {
                        let factor = AICG_CONTACT_AVE / mean_c;
                        for c in &mut raw_contacts {
                            c.3 *= factor;
                        }
                    }
                }
            }
        }
    }

    let mut pairs = Vec::with_capacity(raw_contacts.len());
    for (i, j, r0, e_raw) in raw_contacts {
        if idr_nonlocal.contains(&i) || idr_nonlocal.contains(&j) {
            continue;
        }
        let eps = match options.model {
            ProteinModel::Aicg2Plus => e_raw * CAL_TO_KJ,
            ProteinModel::ClementiGo => {
                CCGO_NATIVE_EPSILON * options.ccgo_contact_scale * CAL_TO_KJ
            }
        };
        pairs.push(ProteinPair {
            i,
            j,
            func: 2,
            eq_nm: r0 * 0.1,
            coef_kj: eps,
        });
    }

    Ok(Protein {
        name: name.to_string(),
        beads,
        chains,
        bonds,
        angles,
        dihedrals,
        pairs,
        model: options.model,
        source_atoms: structure.atoms.len(),
        source_residues: structure.residues.len(),
        warnings,
    })
}

/// Find all native Gō contacts `(i_1based, j_1based, r0_angstrom, raw_energy_kcal)`
/// sorted deterministically by `(i, j)`.
///
/// Uses a 3D cell grid sized to `2 * max_residue_radius + 6.5 Å` so every
/// possible contacting pair lies in the 27-cell neighborhood, making the search
/// strictly exact and O(N) for large complexes.
fn find_native_contacts(
    residues: &[ResidueHeavy],
    model: ProteinModel,
) -> Vec<(usize, usize, f64, f64)> {
    let n = residues.len();
    if n < 2 {
        return Vec::new();
    }

    let max_rad = residues.iter().map(|r| r.radius).fold(0.0_f64, f64::max);
    let cell = (2.0 * max_rad + AICG_GO_ATOMIC_CUTOFF).max(10.0);
    let inv_cell = 1.0 / cell;

    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for r in residues {
        for a in 0..3 {
            lo[a] = lo[a].min(r.center[a]);
            hi[a] = hi[a].max(r.center[a]);
        }
    }
    let dim = |a: usize| (((hi[a] - lo[a]) * inv_cell).floor() as usize + 1).max(1);
    let (nx, ny, nz) = (dim(0), dim(1), dim(2));

    // Cap grid memory for sparse/stretched coordinates.
    let use_grid = nx.saturating_mul(ny).saturating_mul(nz) <= 2_000_000 && n >= 32;

    if !use_grid {
        let per_i: Vec<Vec<(usize, usize, f64, f64)>> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut local = Vec::new();
                let ri = &residues[i];
                for (j, rj) in residues.iter().enumerate().skip(i + 1) {
                    if ri.chain_idx == rj.chain_idx && j < i + 4 {
                        continue;
                    }
                    if is_native_contact(ri, rj) {
                        let r0 = dist(ri.ca_pos, rj.ca_pos);
                        let e = if model == ProteinModel::Aicg2Plus {
                            aicg_raw_energy(&ri.atoms, &rj.atoms)
                        } else {
                            0.0
                        };
                        local.push((i + 1, j + 1, r0, e));
                    }
                }
                local
            })
            .collect();
        return per_i.into_iter().flatten().collect();
    }

    let ncells = nx * ny * nz;
    let mut cell_of = vec![[0usize; 3]; n];
    let mut counts = vec![0u32; ncells];
    for (i, r) in residues.iter().enumerate() {
        let cx = (((r.center[0] - lo[0]) * inv_cell) as usize).min(nx - 1);
        let cy = (((r.center[1] - lo[1]) * inv_cell) as usize).min(ny - 1);
        let cz = (((r.center[2] - lo[2]) * inv_cell) as usize).min(nz - 1);
        cell_of[i] = [cx, cy, cz];
        counts[(cz * ny + cy) * nx + cx] += 1;
    }
    let mut offsets = vec![0u32; ncells + 1];
    let mut acc = 0u32;
    for c in 0..ncells {
        offsets[c] = acc;
        acc += counts[c];
        counts[c] = offsets[c];
    }
    offsets[ncells] = acc;
    let mut sorted = vec![0u32; n];
    for (i, &[cx, cy, cz]) in cell_of.iter().enumerate() {
        let c = (cz * ny + cy) * nx + cx;
        let slot = counts[c] as usize;
        counts[c] += 1;
        sorted[slot] = i as u32;
    }

    let per_i: Vec<Vec<(usize, usize, f64, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut local = Vec::new();
            let ri = &residues[i];
            let [cx, cy, cz] = cell_of[i];
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
                    let start = offsets[row + x0] as usize;
                    let end = offsets[row + x1 + 1] as usize;
                    for &j_u32 in &sorted[start..end] {
                        let j = j_u32 as usize;
                        if j <= i {
                            continue;
                        }
                        let rj = &residues[j];
                        if ri.chain_idx == rj.chain_idx && j < i + 4 {
                            continue;
                        }
                        if is_native_contact(ri, rj) {
                            let r0 = dist(ri.ca_pos, rj.ca_pos);
                            let e = if model == ProteinModel::Aicg2Plus {
                                aicg_raw_energy(&ri.atoms, &rj.atoms)
                            } else {
                                0.0
                            };
                            local.push((i + 1, j + 1, r0, e));
                        }
                    }
                }
            }
            local.sort_unstable_by_key(|c| c.1);
            local
        })
        .collect();

    per_i.into_iter().flatten().collect()
}

/// Render a [`Protein`] as a GENESIS/GROMACS `.itp` `[ moleculetype ]` file,
/// matching `genesis_cg_tool/src/lib/parser_top.jl` column formatting.
pub fn render_protein_itp(protein: &Protein, moleculetype_name: &str) -> String {
    let est = protein.beads.len() * 64
        + protein.bonds.len() * 64
        + protein.angles.len() * 82
        + protein.dihedrals.len() * 92
        + protein.pairs.len() * 64
        + 512;
    let mut s = String::with_capacity(est);

    // [ moleculetype ]
    s.push_str("[ moleculetype ]\n");
    s.push_str(";name            nrexcl\n");
    let _ = writeln!(s, "{moleculetype_name:<16} {:>6}", 3);
    s.push('\n');

    // [ atoms ]
    s.push_str("[ atoms ]\n");
    s.push_str(";       nr type     resnr  res atom   cg   charge     mass\n");
    let mut cur_chain = usize::MAX;
    for bead in &protein.beads {
        if bead.chain_idx != cur_chain {
            cur_chain = bead.chain_idx;
            let _ = writeln!(
                s,
                "; +INFO+ CHAIN: {:>6}     SEGNAME: {:>5}",
                bead.chain_idx + 1,
                bead.seg_name
            );
        }
        let _ = writeln!(
            s,
            "{:>10}{:>5}{:>10}{:>5}{:>5}{:>5} {:>8.3} {:>8.3}",
            bead.atom_id,
            bead.atom_type,
            bead.res_id,
            bead.residue,
            bead.atom,
            1,
            bead.charge,
            bead.mass
        );
    }
    s.push('\n');

    // [ bonds ]
    if !protein.bonds.is_empty() {
        s.push_str("[ bonds ]\n");
        s.push_str(";        i         j    f                eq              coef\n");
        for b in &protein.bonds {
            let _ = writeln!(
                s,
                "{:>10}{:>10}{:>5}{}{}",
                b.i,
                b.j,
                b.func,
                fmt_sci(b.eq_nm, 18),
                fmt_sci(b.coef_kj, 18)
            );
        }
        s.push('\n');
    }

    // [ angles ]
    if !protein.angles.is_empty() {
        s.push_str("[ angles ]\n");
        s.push_str(
            ";        i         j         k    f             eq           coef              w\n",
        );
        for a in &protein.angles {
            match a.func {
                1 => {
                    let _ = writeln!(
                        s,
                        "{:>10}{:>10}{:>10}{:>5}{}{}",
                        a.i,
                        a.j,
                        a.k,
                        a.func,
                        fmt_sci(a.eq, 15),
                        fmt_sci(a.coef_kj, 15)
                    );
                }
                21 => {
                    let _ = writeln!(
                        s,
                        "{:>10}{:>10}{:>10}{:>5}{}{}{}",
                        a.i,
                        a.j,
                        a.k,
                        a.func,
                        fmt_sci(a.eq, 15),
                        fmt_sci(a.coef_kj, 15),
                        fmt_sci(a.width, 15)
                    );
                }
                22 => {
                    let _ = writeln!(s, "{:>10}{:>10}{:>10}{:>5}", a.i, a.j, a.k, a.func);
                }
                _ => {}
            }
        }
        s.push('\n');
    }

    // [ dihedrals ]
    if !protein.dihedrals.is_empty() {
        s.push_str("[ dihedrals ]\n");
        s.push_str(";        i         j         k         l    f             eq           coef            w/n\n");
        for d in &protein.dihedrals {
            match d.func {
                1 | 32 | 33 => {
                    let _ = writeln!(
                        s,
                        "{:>10}{:>10}{:>10}{:>10}{:>5}{}{}{:>15}",
                        d.i,
                        d.j,
                        d.k,
                        d.l,
                        d.func,
                        fmt_sci(d.eq_deg, 15),
                        fmt_sci(d.coef_kj, 15),
                        d.multiplicity
                    );
                }
                21 | 41 | 43 => {
                    let _ = writeln!(
                        s,
                        "{:>10}{:>10}{:>10}{:>10}{:>5}{}{}{}",
                        d.i,
                        d.j,
                        d.k,
                        d.l,
                        d.func,
                        fmt_sci(d.eq_deg, 15),
                        fmt_sci(d.coef_kj, 15),
                        fmt_sci(d.width, 15)
                    );
                }
                22 | 52 => {
                    let _ = writeln!(
                        s,
                        "{:>10}{:>10}{:>10}{:>10}{:>5}",
                        d.i, d.j, d.k, d.l, d.func
                    );
                }
                _ => {}
            }
        }
        s.push('\n');
    }

    // [ pairs ]
    if !protein.pairs.is_empty() {
        s.push_str("[ pairs ]\n");
        s.push_str(";        i         j         f             eq           coef\n");
        for p in &protein.pairs {
            let _ = writeln!(
                s,
                "{:>10}{:>10}{:>10}{}{}",
                p.i,
                p.j,
                p.func,
                fmt_sci(p.eq_nm, 15),
                fmt_sci(p.coef_kj, 15)
            );
        }
        s.push('\n');
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::structure::parse_pdb;

    /// Build a synthetic 10-residue ideal α-helix (all-atom N, CA, C, O, CB) for
    /// self-contained AICG2+ and Clementi Gō tests.
    fn synthetic_helix_pdb() -> String {
        let seq = [
            "ALA", "LEU", "LYS", "GLU", "ALA", "VAL", "ASP", "ARG", "PHE", "GLY",
        ];
        let mut s = String::new();
        let mut serial = 1;
        for (i, &res) in seq.iter().enumerate() {
            let idx = i + 1;
            let theta = i as f64 * (100.0_f64).to_radians();
            let r = 2.3;
            let z = i as f64 * 1.5;
            let ca = [r * theta.cos(), r * theta.sin(), z];
            let n = [
                (r - 0.3) * (theta - 0.3).cos(),
                (r - 0.3) * (theta - 0.3).sin(),
                z - 0.6,
            ];
            let c = [
                (r - 0.2) * (theta + 0.3).cos(),
                (r - 0.2) * (theta + 0.3).sin(),
                z + 0.6,
            ];
            let o = [
                (r + 0.5) * (theta + 0.3).cos(),
                (r + 0.5) * (theta + 0.3).sin(),
                z + 1.4,
            ];
            let mut push = |aname: &str, p: [f64; 3]| {
                let _ = writeln!(
                    s,
                    "ATOM  {:>5} {:<4} {:<3} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00",
                    serial, aname, res, idx, p[0], p[1], p[2]
                );
                serial += 1;
            };
            push("N", n);
            push("CA", ca);
            push("C", c);
            push("O", o);
            if res != "GLY" {
                let cb = [(r + 1.5) * theta.cos(), (r + 1.5) * theta.sin(), z];
                push("CB", cb);
            }
        }
        s.push_str("TER\nEND\n");
        s
    }

    #[test]
    fn selection_and_charge_parsers() {
        assert_eq!(
            parse_selection("1-3, 5 to 7, 10").unwrap(),
            vec![1, 2, 3, 5, 6, 7, 10]
        );
        assert!(parse_selection("5-2").is_err());
        assert!(parse_selection("0").is_err());

        let charges = parse_charge_table("# respac\n1 0.5\n3 -0.75\n").unwrap();
        assert_eq!(charges, vec![(1, 0.5), (3, -0.75)]);
    }

    #[test]
    fn builds_aicg2p_and_clementi_on_synthetic_helix() {
        let structure = parse_pdb(&synthetic_helix_pdb()).unwrap();
        let p = build_protein("helix", &structure, &ProteinOptions::default()).unwrap();
        assert_eq!(p.beads.len(), 10);
        assert_eq!(p.chains.len(), 1);
        assert_eq!(p.bonds.len(), 9);
        // 8 Gaussian 1-3 + 8 flexible-local angles.
        assert_eq!(p.angles.len(), 16);
        // 7 Gaussian 1-4 + 7 flexible-local dihedrals.
        assert_eq!(p.dihedrals.len(), 14);
        // Helical i..i+4 contacts exist (< 6.5 Å).
        assert!(!p.pairs.is_empty());
        // Net charge: LYS(+1) + GLU(-1) + ASP(-1) + ARG(+1) = 0.0.
        assert!(p.net_charge().abs() < 1e-9);

        let itp = p.to_itp();
        assert!(itp.contains("[ moleculetype ]"));
        assert!(itp.contains("9.2383E+04"));
        assert!(itp.contains("   21 "));
        assert!(itp.contains("   22\n"));
        assert!(itp.contains("   41 "));
        assert!(itp.contains("   52\n"));
        assert!(itp.contains("[ pairs ]"));

        // IDR removes local Gaussian terms for beads 1..=4 and pairs touching 1..=4.
        let idr_opts = ProteinOptions {
            idr_local: vec![1, 2, 3, 4],
            idr_nonlocal: vec![1, 2, 3, 4],
            charges: vec![(1, 0.25)],
            ..ProteinOptions::default()
        };
        let p_idr = build_protein("helix_idr", &structure, &idr_opts).unwrap();
        assert!((p_idr.beads[0].charge - 0.25).abs() < 1e-9);
        assert!(
            p_idr
                .angles
                .iter()
                .filter(|a| a.func == 21)
                .all(|a| a.i > 4)
        );
        assert!(p_idr.pairs.iter().all(|c| c.i > 4 && c.j > 4));

        // Clementi Gō model.
        let ccgo = build_protein(
            "helix_go",
            &structure,
            &ProteinOptions {
                model: ProteinModel::ClementiGo,
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        assert_eq!(ccgo.angles.len(), 8);
        assert_eq!(ccgo.dihedrals.len(), 14); // 7 (n=1) + 7 (n=3)
        assert!(ccgo.to_itp().contains("8.3680E+04"));
    }
}
