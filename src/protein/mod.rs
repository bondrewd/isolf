//! Coarse-grained protein representation (AICG2+ / Clementi Gō) for GENESIS RESIDCG.
//!
//! Reads all-atom PDB or mmCIF files ([`structure`]), builds the exact AICG2+
//! or Clementi Gō topology ([`aicg2p`]) using verbatim force-field tables
//! ([`params`]), and renders standalone coordinate and topology outputs or
//! supplies templates for membrane, vesicle, and tubule embedding
//! ([`crate::embed`]).

pub mod aicg2p;
pub mod params;
pub mod structure;

use std::fmt::Write;
use std::ops::Range;
use std::path::Path;

use crate::error::BuildError;
use crate::membrane::Particle;
use crate::protein::aicg2p::{ProteinModel, ProteinOptions, build_protein, render_protein_itp};
use crate::protein::structure::{AllAtomStructure, read_structure};
use crate::voronoi::{cross3, dot3, norm3, rotate_with_sin_cos};

/// One coarse-grained Cα bead of a protein.
#[derive(Debug, Clone, PartialEq)]
pub struct ProteinBead {
    /// Canonical 3-letter amino-acid name (e.g. `"ALA"`).
    pub residue: String,
    /// Atom name (`"CA"`).
    pub atom: String,
    /// Atom type in `aicg2p.itp` (e.g. `"ALA"`).
    pub atom_type: String,
    /// 1-based residue number.
    pub res_id: usize,
    /// 1-based atom number.
    pub atom_id: usize,
    /// 0-based chain index within the protein.
    pub chain_idx: usize,
    /// PDB chain identifier.
    pub chain_id: char,
    /// Segment identifier (from PDB, if any).
    pub seg_name: String,
    /// Cartesian position (nm).
    pub position: [f64; 3],
    /// Charge (e).
    pub charge: f64,
    /// Mass (g/mol).
    pub mass: f64,
    /// Excluded-volume radius (`0.5 * rmin`, nm).
    pub radius: f64,
}

/// One polypeptide chain within a [`Protein`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinChain {
    /// PDB chain identifier.
    pub chain_id: char,
    /// Segment identifier.
    pub seg_name: String,
    /// Contiguous bead index range into [`Protein::beads`].
    pub range: Range<usize>,
}

/// A bonded Cα–Cα pair in `[ bonds ]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinBond {
    pub i: usize,
    pub j: usize,
    pub func: u8,
    pub eq_nm: f64,
    pub coef_kj: f64,
}

/// A 1–3 angle term in `[ angles ]` (`func = 21` Gaussian, `22` flexible, `1` harmonic).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinAngle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub func: u8,
    pub eq: f64,
    pub coef_kj: f64,
    pub width: f64,
}

/// A 1–4 dihedral term in `[ dihedrals ]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinDihedral {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub func: u8,
    pub eq_deg: f64,
    pub coef_kj: f64,
    pub width: f64,
    pub multiplicity: u8,
}

/// A Gō native contact in `[ pairs ]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinPair {
    pub i: usize,
    pub j: usize,
    pub func: u8,
    pub eq_nm: f64,
    pub coef_kj: f64,
}

/// A coarse-grained protein with coordinates (nm) and full GENESIS topology.
#[derive(Debug, Clone, PartialEq)]
pub struct Protein {
    /// Molecule/topology name (used in `[ moleculetype ]` and output filenames).
    pub name: String,
    /// Coarse-grained Cα beads.
    pub beads: Vec<ProteinBead>,
    /// Polypeptide chain ranges.
    pub chains: Vec<ProteinChain>,
    /// Covalent backbone bonds.
    pub bonds: Vec<ProteinBond>,
    /// Local angle potentials.
    pub angles: Vec<ProteinAngle>,
    /// Local dihedral potentials.
    pub dihedrals: Vec<ProteinDihedral>,
    /// Native Gō contact pairs.
    pub pairs: Vec<ProteinPair>,
    /// Force-field model used to build this protein.
    pub model: ProteinModel,
    /// Number of all-atom entries read from the input structure.
    pub source_atoms: usize,
    /// Number of protein residues read from the input structure.
    pub source_residues: usize,
    /// Non-fatal build warnings (e.g. abnormal CA–CA bond lengths).
    pub warnings: Vec<String>,
}

impl Protein {
    /// Build a coarse-grained protein from an already-parsed [`AllAtomStructure`].
    ///
    /// # Errors
    ///
    /// Propagates [`BuildError`] from [`build_protein`].
    pub fn from_structure(
        name: &str,
        structure: &AllAtomStructure,
        options: &ProteinOptions,
    ) -> Result<Self, BuildError> {
        build_protein(name, structure, options)
    }

    /// Read a PDB or mmCIF file and build a coarse-grained protein.
    ///
    /// If `name_override` is `None`, the molecule name defaults to `<stem>_cg`
    /// sanitized to alphanumeric/underscore characters.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError`] if the file cannot be read/parsed or contains no
    /// coarse-grainable protein residues.
    pub fn from_file(
        path: &Path,
        name_override: Option<&str>,
        options: &ProteinOptions,
    ) -> Result<Self, BuildError> {
        let structure = read_structure(path)?;
        let name = match name_override {
            Some(n) => sanitize_mol_name(n),
            None => {
                let stem = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("protein");
                sanitize_mol_name(&format!("{stem}_cg"))
            }
        };
        Self::from_structure(&name, &structure, options)
    }

    /// Geometric centroid of all beads (nm).
    pub fn centroid(&self) -> [f64; 3] {
        if self.beads.is_empty() {
            return [0.0; 3];
        }
        let inv = 1.0 / self.beads.len() as f64;
        let mut sum = [0.0; 3];
        for b in &self.beads {
            sum[0] += b.position[0];
            sum[1] += b.position[1];
            sum[2] += b.position[2];
        }
        [sum[0] * inv, sum[1] * inv, sum[2] * inv]
    }

    /// Radius of gyration across all beads (nm).
    pub fn radius_of_gyration(&self) -> f64 {
        range_rg(&self.beads, 0..self.beads.len())
    }

    /// Radius of gyration of one polypeptide chain (nm).
    pub fn chain_rg(&self, chain_idx: usize) -> f64 {
        self.chains
            .get(chain_idx)
            .map_or(0.0, |c| range_rg(&self.beads, c.range.clone()))
    }

    /// Radius of the circumsphere about the centroid (nm).
    pub fn circumsphere_radius(&self) -> f64 {
        let c = self.centroid();
        self.beads
            .iter()
            .map(|b| {
                let dx = b.position[0] - c[0];
                let dy = b.position[1] - c[1];
                let dz = b.position[2] - c[2];
                dx * dx + dy * dy + dz * dz
            })
            .fold(0.0_f64, f64::max)
            .sqrt()
    }

    /// Axis-aligned bounding box `(min, max)` in nm.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for b in &self.beads {
            for a in 0..3 {
                lo[a] = lo[a].min(b.position[a]);
                hi[a] = hi[a].max(b.position[a]);
            }
        }
        if self.beads.is_empty() {
            ([0.0; 3], [0.0; 3])
        } else {
            (lo, hi)
        }
    }

    /// Total net charge across all beads (e).
    pub fn net_charge(&self) -> f64 {
        self.beads.iter().map(|b| b.charge).sum()
    }

    /// Total net charge of one polypeptide chain (e).
    pub fn chain_charge(&self, chain_idx: usize) -> f64 {
        self.chains.get(chain_idx).map_or(0.0, |c| {
            self.beads[c.range.clone()].iter().map(|b| b.charge).sum()
        })
    }

    /// Translate every bead by `shift` (nm) in place.
    pub fn translate(&mut self, shift: [f64; 3]) {
        for b in &mut self.beads {
            b.position[0] += shift[0];
            b.position[1] += shift[1];
            b.position[2] += shift[2];
        }
    }

    /// Rotate every bead about `origin` around `axis` by `angle_rad` (Rodrigues).
    pub fn rotate_axis(&mut self, origin: [f64; 3], axis: [f64; 3], angle_rad: f64) {
        if angle_rad.abs() < 1e-12 {
            return;
        }
        let m = norm3(axis);
        if m < 1e-12 {
            return;
        }
        let u = [axis[0] / m, axis[1] / m, axis[2] / m];
        let (sin_a, cos_a) = angle_rad.sin_cos();
        for b in &mut self.beads {
            let rel = [
                b.position[0] - origin[0],
                b.position[1] - origin[1],
                b.position[2] - origin[2],
            ];
            let r = rotate_with_sin_cos(rel, u, sin_a, cos_a);
            b.position = [origin[0] + r[0], origin[1] + r[1], origin[2] + r[2]];
        }
    }

    /// Map local coordinates `(x, y, z)` into an orthonormal frame `(u, v, n)`
    /// anchored at `origin`: `p' = origin + x·u + y·v + z·n`.
    pub fn apply_frame(&mut self, origin: [f64; 3], u: [f64; 3], v: [f64; 3], n: [f64; 3]) {
        for b in &mut self.beads {
            let [x, y, z] = b.position;
            b.position = [
                origin[0] + x * u[0] + y * v[0] + z * n[0],
                origin[1] + x * u[1] + y * v[1] + z * n[1],
                origin[2] + x * u[2] + y * v[2] + z * n[2],
            ];
        }
    }

    /// Cubic periodic box spanning the protein plus `padding` nm on each side.
    pub fn box_with_padding(&self, padding: f64) -> [f64; 3] {
        let (lo, hi) = self.bounding_box();
        let side = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(hi[2] - lo[2]) + 2.0 * padding;
        [
            side.max(2.0 * padding),
            side.max(2.0 * padding),
            side.max(2.0 * padding),
        ]
    }

    /// Convert beads to [`Particle`] records (1 residue id per bead, matching
    /// GENESIS `.gro` / `.pdb` / `.psf` conventions).
    pub fn to_particles(&self, start_res_id: usize, start_atom_id: usize) -> Vec<Particle> {
        self.beads
            .iter()
            .enumerate()
            .map(|(k, b)| Particle {
                residue: b.residue.clone(),
                bead: b.atom.clone(),
                molecule_id: start_res_id + k,
                atom_id: start_atom_id + k,
                position: b.position,
                charge: b.charge,
                mass: b.mass,
            })
            .collect()
    }

    /// Render this protein's `[ moleculetype ]` `.itp` file.
    pub fn to_itp(&self) -> String {
        render_protein_itp(self, &self.name)
    }

    /// Render a standalone GENESIS `.top` file referencing `./aicg2p.itp` and
    /// `./<name>.itp`.
    pub fn to_top(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "; AICG2+ protein force field parameters");
        let _ = writeln!(s, "#include \"./aicg2p.itp\"");
        let _ = writeln!(s);
        let _ = writeln!(s, "; Protein molecule topology");
        let _ = writeln!(s, "#include \"./{}.itp\"", self.name);
        let _ = writeln!(s);
        let _ = writeln!(s, "[ system ]");
        let _ = writeln!(s, "{}", self.name);
        let _ = writeln!(s);
        let _ = writeln!(s, "[ molecules ]");
        let _ = writeln!(s, "{} 1", self.name);
        let _ = writeln!(s);
        let _ = writeln!(s, "[ cg_ele_chain_pairs ]");
        let _ = writeln!(s, "ON 1 - 1 : 1 - 1");
        s
    }

    /// Render a standalone GROMACS `.gro` file with the given `box_size` (nm).
    pub fn to_gro(&self, box_size: [f64; 3]) -> String {
        let mut s = String::with_capacity(self.beads.len() * 72 + 128);
        let _ = writeln!(s, "CG model {}, t = {:>16.3} ", self.name, 0.0);
        let _ = writeln!(s, "{:>12} ", self.beads.len());
        for b in &self.beads {
            let _ = writeln!(
                s,
                "{:>5}{:>5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}  0.0000  0.0000  0.0000 ",
                b.res_id % 100000,
                b.residue,
                b.atom,
                b.atom_id % 100000,
                b.position[0],
                b.position[1],
                b.position[2],
            );
        }
        let _ = writeln!(
            s,
            "{:>15.4}{:>15.4}{:>15.4} ",
            box_size[0], box_size[1], box_size[2]
        );
        s.push('\n');
        s
    }

    /// Render a chain-aware CHARMM `.psf` file for this protein alone.
    pub fn to_psf(&self) -> String {
        let particles = self.to_particles(1, 1);
        let chains: Vec<(String, Range<usize>)> = self
            .chains
            .iter()
            .enumerate()
            .map(|(i, c)| (protein_segname(i, &c.seg_name), c.range.clone()))
            .collect();
        crate::psf::render_chains(&particles, &chains, &self.name)
    }

    /// Render a chain-aware `.pdb` file (coordinates in Å, `TER` between chains).
    pub fn to_pdb(&self) -> String {
        let mut s = String::with_capacity(self.beads.len() * 82 + 128);
        let _ = writeln!(s, "TITLE     {:<70}", self.name);
        for ch in &self.chains {
            let cid = if ch.chain_id == ' ' { 'A' } else { ch.chain_id };
            for b in &self.beads[ch.range.clone()] {
                let _ = writeln!(
                    s,
                    "ATOM  {:>5} {:>4} {:<4}{}{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00               ",
                    b.atom_id % 100000,
                    b.atom,
                    b.residue,
                    cid,
                    b.res_id % 10000,
                    b.position[0] * 10.0,
                    b.position[1] * 10.0,
                    b.position[2] * 10.0,
                );
            }
            s.push_str("TER\n");
        }
        s.push_str("END\n");
        s
    }

    /// Render a CHARMM `.crd` file for this protein alone.
    pub fn to_crd(&self) -> String {
        let particles = self.to_particles(1, 1);
        crate::crd::render(&particles, &self.name, self.beads.len())
    }

    /// Render an mmCIF `.cif` file for this protein alone.
    pub fn to_cif(&self) -> String {
        let particles = self.to_particles(1, 1);
        crate::cif::render(&particles, &self.name)
    }
}

/// Format a 4-character CHARMM PSF segment name for protein chain `idx`.
pub(crate) fn protein_segname(chain_idx: usize, existing: &str) -> String {
    let trimmed = existing.trim();
    if !trimmed.is_empty() && trimmed.len() <= 4 {
        return trimmed.to_string();
    }
    if chain_idx < 26 {
        format!("PRO{}", (b'A' + chain_idx as u8) as char)
    } else {
        format!("P{:03}", (chain_idx + 1) % 1000)
    }
}

/// Sanitize a molecule name for GROMACS/GENESIS `[ moleculetype ]` (no whitespace).
pub fn sanitize_mol_name(raw: &str) -> String {
    let clean: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if clean.is_empty() {
        "protein_cg".to_string()
    } else {
        clean
    }
}

fn range_rg(beads: &[ProteinBead], range: Range<usize>) -> f64 {
    let slice = &beads[range];
    if slice.is_empty() {
        return 0.0;
    }
    let inv = 1.0 / slice.len() as f64;
    let mut c = [0.0; 3];
    for b in slice {
        c[0] += b.position[0];
        c[1] += b.position[1];
        c[2] += b.position[2];
    }
    c[0] *= inv;
    c[1] *= inv;
    c[2] *= inv;
    let sq_sum: f64 = slice
        .iter()
        .map(|b| {
            let dx = b.position[0] - c[0];
            let dy = b.position[1] - c[1];
            let dz = b.position[2] - c[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sq_sum * inv).sqrt()
}

/// Build an orthonormal tangent frame `(u, v, n)` from a unit normal `n`.
pub(crate) fn orthonormal_frame(normal: [f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let m = norm3(normal).max(1e-12);
    let n = [normal[0] / m, normal[1] / m, normal[2] / m];
    // Pick a reference vector not parallel to `n` (prefer +x near +z pole so a
    // protein placed at the north pole keeps its exact PDB x/y orientation).
    let ref_vec = if n[2].abs() > 0.99 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let v_raw = cross3(n, ref_vec);
    let vm = norm3(v_raw).max(1e-12);
    let v = [v_raw[0] / vm, v_raw[1] / vm, v_raw[2] / vm];
    let u = cross3(v, n);
    // Ensure right-handedness: (u × v) · n > 0.
    if dot3(cross3(u, v), n) < 0.0 {
        ([-u[0], -u[1], -u[2]], v, n)
    } else {
        (u, v, n)
    }
}
