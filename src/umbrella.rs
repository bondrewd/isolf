//! Umbrella sampling setup, lipid–protein interaction `.itp` generation,
//! GENESIS CV log extraction, and MBAR/BO campaign manifest utilities.
//!
//! Designed for parameterizing iSoLFv2 lipid–protein interactions via amino-acid
//! solvation PMFs across a flat lipid bilayer centered at `z = 0`:
//!
//! - **Single-bead probe placement** without carving lipids, so every window
//!   shares the exact same lipid count and Hamiltonian (required by MBAR/WHAM).
//! - **Decoupled `lipid_protein.itp`** carrying `[ cg_ISOLF_parameters ]`
//!   (WCA + cubic polynomial switch) plus optional `[ cg_WCA_parameters ]` /
//!   `[ cg_LJ_parameters ]` head-group interactions, so Bayesian optimization
//!   iterations only rewrite `lipid_protein.itp` without rebuilding windows.
//! - **GENESIS `atdyn` control files** using signed COM z-distance restraints
//!   (`function1 = DISTMASS`, `direction1 = Z2`, `U = k * (dz - z0)^2`) in
//!   either independent-window mode or Replica-Exchange Umbrella Sampling
//!   (`REUS`), paired with `mbar_analysis` (`mbar.inp`) and a fast `RESTR_CVS1`
//!   log extractor.

use std::fmt::Write as _;

use crate::error::BuildError;
use crate::force_field::{ForceField, Interaction, LJ_CUTOFF_RATIO, POTENTIAL_MINIMUM_FACTOR};
use crate::membrane::{Membrane, Particle, molecule_ranges};
use crate::protein::{Protein, params};

/// Default iSoLF well depth `epsilon` (kJ/mol) for an unparameterized amino-acid
/// probe before optimization.
pub const DEFAULT_PROBE_EPSILON: f64 = 1.5;

/// Default iSoLF polynomial-switch attraction width `omega` (nm).
pub const DEFAULT_PROBE_OMEGA: f64 = 1.0;

/// Default harmonic umbrella force constant `k` in kcal·mol⁻¹·Å⁻²
/// (`U = k * (dz - z0)^2`, matching GENESIS `atdyn` & `mbar_analysis`).
pub const DEFAULT_UMBRELLA_K: f64 = 0.5;

/// Default z-box length (nm) for bilayer umbrella sampling (`200 Å`, matching
/// Ugarte La Torre & Takada 2021).
pub const DEFAULT_UMBRELLA_BOX_Z: f64 = 20.0;

/// One 1-body entry in `[ cg_ISOLF_parameters ]` (WCA repulsion + cubic
/// polynomial switch attraction combined via Lorentz–Berthelot / arithmetic
/// rules with all other iSoLF beads).
#[derive(Debug, Clone, PartialEq)]
pub struct IsolfBeadParam {
    /// Atom-type name (e.g. `"TRP"`).
    pub name: String,
    /// Well depth `epsilon` (kJ/mol, must be `>= 0`).
    pub epsilon: f64,
    /// Repulsive diameter `sigma` (nm, must be `> 0`).
    pub sigma: f64,
    /// Polynomial-switch attraction width `omega` (nm, must be `> 0` when `epsilon > 0`).
    pub omega: f64,
}

/// One explicit pair entry in `[ cg_WCA_parameters ]` (purely repulsive WCA).
#[derive(Debug, Clone, PartialEq)]
pub struct WcaPairParam {
    pub first: String,
    pub second: String,
    /// Well depth `epsilon` (kJ/mol, `>= 0`).
    pub epsilon: f64,
    /// Repulsive diameter `sigma` (nm, `> 0`).
    pub sigma: f64,
}

/// One explicit pair entry in `[ cg_LJ_parameters ]` (attractive LJ with cutoff).
#[derive(Debug, Clone, PartialEq)]
pub struct LjPairParam {
    pub first: String,
    pub second: String,
    /// Well depth `epsilon` (kJ/mol, `>= 0`).
    pub epsilon: f64,
    /// Diameter `sigma` (nm, `> 0`).
    pub sigma: f64,
    /// Cutoff distance (nm, defaults to `2.5 * sigma`).
    pub cutoff: f64,
}

/// Tunable lipid–protein non-bonded interaction parameters rendered into a
/// standalone `lipid_protein.itp` file for GENESIS `RESIDCG`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LipidProteinParams {
    /// Per-atomtype `[ cg_ISOLF_parameters ]` entries (WCA + polynomial switch).
    pub isolf: Vec<IsolfBeadParam>,
    /// Explicit `[ cg_WCA_parameters ]` pairs (e.g. amino acid ↔ lipid head).
    pub wca: Vec<WcaPairParam>,
    /// Explicit `[ cg_LJ_parameters ]` pairs (optional polar/head attraction).
    pub lj: Vec<LjPairParam>,
}

/// Inline CLI overrides for a single probe's `[ cg_ISOLF_parameters ]` entry.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ProbeParamOverrides {
    pub epsilon: Option<f64>,
    pub sigma: Option<f64>,
    pub omega: Option<f64>,
    pub auto_head_wca: bool,
}

impl LipidProteinParams {
    /// Default repulsive diameter `sigma` (nm) for a canonical amino acid derived
    /// from its AICG2+ `rmin` (`sigma = rmin / 2^(1/6)`).
    pub fn default_sigma(aa: &str) -> f64 {
        params::aa_rmin(aa) / POTENTIAL_MINIMUM_FACTOR
    }

    /// Parse a lipid–protein parameter file with `[cg_ISOLF_parameters]` (or
    /// `[isolf]`), `[cg_WCA_parameters]` (or `[wca]`), and `[cg_LJ_parameters]`
    /// (or `[lj]`) sections.
    ///
    /// Lines starting with `#` or `;` are comments. Both whitespace- and
    /// comma-delimited columns are accepted.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::InvalidSelection`] on syntax errors or invalid
    /// numeric ranges.
    pub fn parse(content: &str) -> Result<Self, BuildError> {
        #[derive(Clone, Copy)]
        enum Section {
            None,
            Isolf,
            Wca,
            Lj,
        }
        let mut out = Self::default();
        let mut section = Section::None;

        for (line_no, raw) in content.lines().enumerate() {
            let line = raw.split([';', '#']).next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                let hdr = line[1..line.len() - 1].trim().to_ascii_lowercase();
                section = match hdr.as_str() {
                    "cg_isolf_parameters" | "isolf" => Section::Isolf,
                    "cg_wca_parameters" | "wca" => Section::Wca,
                    "cg_lj_parameters" | "lj" => Section::Lj,
                    other => {
                        return Err(BuildError::InvalidSelection(format!(
                            "params line {}: unknown section '[{other}]' \
                             (expected [cg_ISOLF_parameters], [cg_WCA_parameters], or [cg_LJ_parameters])",
                            line_no + 1
                        )));
                    }
                };
                continue;
            }

            let toks: Vec<&str> = line
                .split(|c: char| c.is_whitespace() || c == ',')
                .filter(|t| !t.is_empty())
                .collect();
            let err = |msg: &str| {
                BuildError::InvalidSelection(format!("params line {}: {msg}", line_no + 1))
            };
            let parse_num = |s: &str, what: &str| -> Result<f64, BuildError> {
                s.parse::<f64>()
                    .ok()
                    .filter(|v| v.is_finite())
                    .ok_or_else(|| err(&format!("invalid {what} '{s}'")))
            };

            match section {
                Section::None => {
                    return Err(err(
                        "parameter entry outside any section header (add [cg_ISOLF_parameters], [cg_WCA_parameters], or [cg_LJ_parameters])",
                    ));
                }
                Section::Isolf => {
                    if toks.len() != 4 {
                        return Err(err("expected 4 columns: <name> <epsilon> <sigma> <omega>"));
                    }
                    let name = toks[0].to_ascii_uppercase();
                    let epsilon = parse_num(toks[1], "epsilon")?;
                    let sigma = parse_num(toks[2], "sigma")?;
                    let omega = parse_num(toks[3], "omega")?;
                    validate_isolf_entry(&name, epsilon, sigma, omega)?;
                    if let Some(existing) = out.isolf.iter_mut().find(|e| e.name == name) {
                        existing.epsilon = epsilon;
                        existing.sigma = sigma;
                        existing.omega = omega;
                    } else {
                        out.isolf.push(IsolfBeadParam {
                            name,
                            epsilon,
                            sigma,
                            omega,
                        });
                    }
                }
                Section::Wca => {
                    if toks.len() != 4 {
                        return Err(err("expected 4 columns: <name1> <name2> <epsilon> <sigma>"));
                    }
                    let first = toks[0].to_ascii_uppercase();
                    let second = toks[1].to_ascii_uppercase();
                    let epsilon = parse_num(toks[2], "epsilon")?;
                    let sigma = parse_num(toks[3], "sigma")?;
                    if epsilon < 0.0 || sigma <= 0.0 {
                        return Err(err("WCA requires epsilon >= 0 and sigma > 0"));
                    }
                    out.upsert_wca(WcaPairParam {
                        first,
                        second,
                        epsilon,
                        sigma,
                    });
                }
                Section::Lj => {
                    if toks.len() != 4 && toks.len() != 5 {
                        return Err(err(
                            "expected 4 or 5 columns: <name1> <name2> <epsilon> <sigma> [cutoff]",
                        ));
                    }
                    let first = toks[0].to_ascii_uppercase();
                    let second = toks[1].to_ascii_uppercase();
                    let epsilon = parse_num(toks[2], "epsilon")?;
                    let sigma = parse_num(toks[3], "sigma")?;
                    let cutoff = if toks.len() == 5 {
                        parse_num(toks[4], "cutoff")?
                    } else {
                        sigma * LJ_CUTOFF_RATIO
                    };
                    if epsilon < 0.0 || sigma <= 0.0 || cutoff <= 0.0 {
                        return Err(err("LJ requires epsilon >= 0, sigma > 0, and cutoff > 0"));
                    }
                    out.upsert_lj(LjPairParam {
                        first,
                        second,
                        epsilon,
                        sigma,
                        cutoff,
                    });
                }
            }
        }
        Ok(out)
    }

    /// Build a complete parameter table for `probe_aas`, starting from an optional
    /// user parameter file, applying inline CLI overrides, and optionally filling
    /// missing amino-acid ↔ lipid-head WCA excluded-volume pairs using
    /// Lorentz–Berthelot combining rules.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::InvalidSelection`] if any override value is out of range.
    pub fn resolve(
        force_field: &ForceField,
        probe_aas: &[&str],
        base: Option<&Self>,
        overrides: ProbeParamOverrides,
    ) -> Result<Self, BuildError> {
        let mut out = base.cloned().unwrap_or_default();

        // Ensure every requested probe amino acid has an entry in `[ cg_ISOLF_parameters ]`
        // (unless the user explicitly omitted it while supplying WCA/LJ pairs for it).
        for &aa in probe_aas {
            let has_explicit_any = out.isolf.iter().any(|e| e.name == aa)
                || out.wca.iter().any(|p| p.first == aa || p.second == aa)
                || out.lj.iter().any(|p| p.first == aa || p.second == aa);

            if let Some(entry) = out.isolf.iter_mut().find(|e| e.name == aa) {
                if let Some(e) = overrides.epsilon {
                    entry.epsilon = e;
                }
                if let Some(s) = overrides.sigma {
                    entry.sigma = s;
                }
                if let Some(w) = overrides.omega {
                    entry.omega = w;
                }
                validate_isolf_entry(&entry.name, entry.epsilon, entry.sigma, entry.omega)?;
            } else if !has_explicit_any
                || overrides.epsilon.is_some()
                || overrides.sigma.is_some()
                || overrides.omega.is_some()
            {
                let epsilon = overrides.epsilon.unwrap_or(DEFAULT_PROBE_EPSILON);
                let sigma = overrides.sigma.unwrap_or_else(|| Self::default_sigma(aa));
                let omega = overrides.omega.unwrap_or(DEFAULT_PROBE_OMEGA);
                validate_isolf_entry(aa, epsilon, sigma, omega)?;
                out.isolf.push(IsolfBeadParam {
                    name: aa.to_string(),
                    epsilon,
                    sigma,
                    omega,
                });
            }
        }

        // Fill Lorentz-Berthelot WCA repulsion between each probe AA and every
        // non-iSoLF lipid head bead (PHA, CHO, ETH, GLC, SRI, PHO, MID) unless
        // that pair already has an explicit WCA or LJ entry.
        if overrides.auto_head_wca {
            let head_beads: Vec<_> = force_field
                .bead_types
                .iter()
                .filter(|b| !matches!(b.interaction, Interaction::Isolf { .. }))
                .collect();
            for &aa in probe_aas {
                let (aa_eps, aa_sig) = out
                    .isolf
                    .iter()
                    .find(|e| e.name == aa)
                    .map(|e| (e.epsilon, e.sigma))
                    .unwrap_or((DEFAULT_PROBE_EPSILON, Self::default_sigma(aa)));
                let wca_eps_aa = if aa_eps > 0.0 {
                    aa_eps
                } else {
                    DEFAULT_PROBE_EPSILON
                };
                for hb in &head_beads {
                    if out.has_pair(aa, &hb.name) {
                        continue;
                    }
                    let eps = (wca_eps_aa * hb.interaction.epsilon()).sqrt();
                    let sig = 0.5 * (aa_sig + hb.interaction.sigma());
                    out.wca.push(WcaPairParam {
                        first: aa.to_string(),
                        second: hb.name.clone(),
                        epsilon: eps,
                        sigma: sig,
                    });
                }
            }
        }

        // Never emit a pair in both WCA and LJ; explicit LJ takes precedence.
        let lj_pairs: Vec<(String, String)> = out
            .lj
            .iter()
            .map(|p| (p.first.clone(), p.second.clone()))
            .collect();
        out.wca.retain(|w| {
            !lj_pairs.iter().any(|(a, b)| {
                (w.first == *a && w.second == *b) || (w.first == *b && w.second == *a)
            })
        });

        Ok(out)
    }

    fn has_pair(&self, a: &str, b: &str) -> bool {
        let eq = |x: &str, y: &str| (x == a && y == b) || (x == b && y == a);
        self.wca.iter().any(|p| eq(&p.first, &p.second))
            || self.lj.iter().any(|p| eq(&p.first, &p.second))
    }

    fn upsert_wca(&mut self, pair: WcaPairParam) {
        if let Some(existing) = self.wca.iter_mut().find(|p| {
            (p.first == pair.first && p.second == pair.second)
                || (p.first == pair.second && p.second == pair.first)
        }) {
            *existing = pair;
        } else {
            self.wca.push(pair);
        }
    }

    fn upsert_lj(&mut self, pair: LjPairParam) {
        if let Some(existing) = self.lj.iter_mut().find(|p| {
            (p.first == pair.first && p.second == pair.second)
                || (p.first == pair.second && p.second == pair.first)
        }) {
            *existing = pair;
        } else {
            self.lj.push(pair);
        }
    }

    /// Render a GENESIS-compatible `lipid_protein.itp` include file.
    pub fn render_itp(&self) -> String {
        let mut s = String::with_capacity(512 + self.isolf.len() * 48 + self.wca.len() * 48);
        s.push_str("; ------------------------------------------------\n");
        s.push_str("; iSoLF lipid-protein interaction parameters\n");
        s.push_str("; ------------------------------------------------\n\n");

        if !self.lj.is_empty() {
            s.push_str("[ cg_LJ_parameters ]\n");
            s.push_str("; name  name  epsilon    sigma  cut-off\n");
            s.push_str(";    -     -   kJ/mol       nm       nm\n");
            for p in &self.lj {
                let _ = writeln!(
                    s,
                    "{:>6}{:>6} {:>8.3} {:>8.3} {:>8.3}",
                    p.first, p.second, p.epsilon, p.sigma, p.cutoff
                );
            }
            s.push('\n');
        }

        if !self.wca.is_empty() {
            s.push_str("[ cg_WCA_parameters ]\n");
            s.push_str("; name  name  epsilon    sigma\n");
            s.push_str(";    -     -   kJ/mol       nm\n");
            for p in &self.wca {
                let _ = writeln!(
                    s,
                    "{:>6}{:>6} {:>8.3} {:>8.3}",
                    p.first, p.second, p.epsilon, p.sigma
                );
            }
            s.push('\n');
        }

        if !self.isolf.is_empty() {
            s.push_str("[ cg_ISOLF_parameters ]\n");
            s.push_str("; name  epsilon    sigma    omega\n");
            s.push_str(";    -   kJ/mol       nm       nm\n");
            for p in &self.isolf {
                let _ = writeln!(
                    s,
                    "{:>6} {:>8.3} {:>8.3} {:>8.3}",
                    p.name, p.epsilon, p.sigma, p.omega
                );
            }
            s.push('\n');
        }
        s
    }
}

fn validate_isolf_entry(
    name: &str,
    epsilon: f64,
    sigma: f64,
    omega: f64,
) -> Result<(), BuildError> {
    if !epsilon.is_finite() || epsilon < 0.0 {
        return Err(BuildError::InvalidSelection(format!(
            "cg_ISOLF_parameters '{name}': epsilon must be >= 0 (got {epsilon})"
        )));
    }
    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(BuildError::InvalidSelection(format!(
            "cg_ISOLF_parameters '{name}': sigma must be > 0 (got {sigma})"
        )));
    }
    if !omega.is_finite() || omega < 0.0 {
        return Err(BuildError::InvalidSelection(format!(
            "cg_ISOLF_parameters '{name}': omega must be >= 0 (got {omega})"
        )));
    }
    if epsilon > 0.0 && omega == 0.0 {
        return Err(BuildError::InvalidSelection(format!(
            "cg_ISOLF_parameters '{name}': omega must be > 0 when epsilon > 0 \
             (use [cg_WCA_parameters] for purely repulsive interactions)"
        )));
    }
    Ok(())
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

/// Center a lipid bilayer's mass-weighted COM at `(0, 0, 0)` (in `z`) and wrap
/// lateral coordinates cleanly around `(0, 0)` in `[-Lx/2, Lx/2) × [-Ly/2, Ly/2)`.
///
/// Setting `z_com = 0` guarantees that a probe placed at `z = z_window` starts
/// at exact signed COM distance `dz = z1 - z2 = z_window` for GENESIS
/// `DISTMASS` + `direction1 = Z2`.
pub fn center_bilayer_at_origin(bilayer: &mut Membrane, min_box_z: f64) {
    if bilayer.particles.is_empty() {
        return;
    }
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

/// Build one umbrella-sampling window coordinate set with the single-bead
/// `probe` at residue 1 / atom 1 (`z = z_nm`) and ALL lipids of `bilayer`
/// preserved intact (residue `2..=N_lipids+1`, atom `2..=N_atoms`).
///
/// The lateral `(x, y)` coordinate of the probe is chosen on a fine 2D grid to
/// maximize clearance from nearby lipid beads in the `z_nm` slab, minimizing
/// initial overlap without deleting any lipids.
pub fn build_window_particles(
    bilayer: &Membrane,
    probe: &Protein,
    z_nm: f64,
    force_field: &ForceField,
) -> Vec<Particle> {
    let bead = &probe.beads[0];
    let [bx, by, _] = bilayer.box_size;
    let (px, py) = best_lateral_xy(bilayer, z_nm, bead.radius, bx, by, force_field);

    let mut out = Vec::with_capacity(1 + bilayer.particles.len());
    out.push(Particle {
        residue: bead.residue.clone(),
        bead: bead.atom.clone(),
        molecule_id: 1,
        atom_id: 1,
        position: [px, py, z_nm],
        charge: bead.charge,
        mass: bead.mass,
    });
    for p in &bilayer.particles {
        let mut lp = p.clone();
        lp.molecule_id += 1;
        lp.atom_id += 1;
        out.push(lp);
    }
    out
}

/// Search a 24×24 lateral grid over `[-bx/2, bx/2) × [-by/2, by/2)` to find the
/// `(x, y)` with maximum minimum clearance from lipid beads near `z_nm`.
fn best_lateral_xy(
    bilayer: &Membrane,
    z_nm: f64,
    probe_r: f64,
    bx: f64,
    by: f64,
    force_field: &ForceField,
) -> (f64, f64) {
    const GRID: usize = 24;
    // Collect lipid beads within 1.5 nm vertically of the target z slice.
    let slab: Vec<([f64; 3], f64)> = bilayer
        .particles
        .iter()
        .filter(|p| (p.position[2] - z_nm).abs() <= 1.5)
        .map(|p| {
            let r = force_field
                .bead_type(&p.bead)
                .map(|bt| 0.5 * POTENTIAL_MINIMUM_FACTOR * bt.interaction.sigma())
                .unwrap_or(0.35);
            (p.position, r + probe_r)
        })
        .collect();
    if slab.is_empty() {
        return (0.0, 0.0);
    }

    let half_x = 0.5 * bx;
    let half_y = 0.5 * by;
    let mut best_xy = (0.0, 0.0);
    let mut best_clearance = f64::NEG_INFINITY;
    for ix in 0..GRID {
        let x = -half_x + ((ix as f64 + 0.5) / GRID as f64) * bx;
        for iy in 0..GRID {
            let y = -half_y + ((iy as f64 + 0.5) / GRID as f64) * by;
            let mut min_clear = f64::INFINITY;
            for &(pos, sum_r) in &slab {
                let mut dx = (x - pos[0]).abs();
                if dx > half_x {
                    dx = bx - dx;
                }
                let mut dy = (y - pos[1]).abs();
                if dy > half_y {
                    dy = by - dy;
                }
                let dz = z_nm - pos[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt() - sum_r;
                if d < min_clear {
                    min_clear = d;
                    if min_clear <= best_clearance {
                        break;
                    }
                }
            }
            if min_clear > best_clearance {
                best_clearance = min_clear;
                best_xy = (x, y);
            }
        }
    }
    best_xy
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

/// Render the shared GENESIS `.top` topology for an umbrella-sampling campaign,
/// including `aicg2p.itp`, `isolf.itp`, `lipid_protein.itp`, and `<probe>.itp`.
pub fn render_umbrella_top(
    probe: &Protein,
    lipid_counts: &[(String, usize)],
    include_prefix: &str,
) -> String {
    let total_mols: usize = 1 + lipid_counts.iter().map(|(_, c)| c).sum::<usize>();
    let mut s = String::new();
    let _ = writeln!(s, "; AICG2+ protein force field parameters");
    let _ = writeln!(s, "#include \"{include_prefix}aicg2p.itp\"");
    let _ = writeln!(s, "; Implicit Solvent Lipid Forcefield (iSoLF)");
    let _ = writeln!(s, "#include \"{include_prefix}isolf.itp\"");
    let _ = writeln!(s, "; Lipid-protein interaction parameters");
    let _ = writeln!(s, "#include \"{include_prefix}lipid_protein.itp\"");
    let _ = writeln!(s, "; Probe amino acid topology");
    let _ = writeln!(s, "#include \"{include_prefix}{}.itp\"", probe.name);
    let _ = writeln!(s);
    let _ = writeln!(s, "[ system ]");
    let _ = writeln!(s, "{} umbrella sampling", probe.name);
    let _ = writeln!(s);
    let _ = writeln!(s, "[ molecules ]");
    let _ = writeln!(s, "{} 1", probe.name);
    for (lipid, count) in lipid_counts {
        let _ = writeln!(s, "{lipid} {count}");
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "[ cg_ele_chain_pairs ]");
    let _ = writeln!(s, "ON 1 - {total_mols} : 1 - {total_mols}");
    s
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
/// Returns [`BuildError::StructureParse`] if any log is missing `RESTR_CVS1`.
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

/// Summary metadata written to `manifest.json` so Bayesian-optimization drivers
/// can discover windows, paths, and force constants without ad-hoc parsing.
#[derive(Debug, Clone)]
pub struct UmbrellaManifest<'a> {
    pub probe: &'a str,
    pub mode: &'a str,
    pub ensemble: &'a str,
    pub temperature_k: f64,
    pub k_kcal_mol_a2: f64,
    pub restraint_function: &'a str,
    pub restraint_direction: &'a str,
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

/// Extract a small set of fields from a previously generated `manifest.json`
/// (used by `--update-params` and `--extract-cv`).
#[derive(Debug, Clone)]
pub struct ParsedManifest {
    pub probe: String,
    pub mode: String,
    pub n_windows: usize,
}

/// Parse the top-level `"probe"`, `"mode"`, and window count from `manifest.json`.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] if required fields are missing.
pub fn parse_manifest(json: &str) -> Result<ParsedManifest, BuildError> {
    let find_str = |key: &str| -> Option<String> {
        let needle = format!("\"{key}\": \"");
        let start = json.find(&needle)? + needle.len();
        let end = json[start..].find('"')?;
        Some(json[start..start + end].to_string())
    };
    let probe = find_str("probe").ok_or_else(|| {
        BuildError::StructureParse("manifest.json: missing 'probe' field".to_string())
    })?;
    let mode = find_str("mode").ok_or_else(|| {
        BuildError::StructureParse("manifest.json: missing 'mode' field".to_string())
    })?;
    let n_windows = json.matches("\"z_angstrom\":").count();
    if n_windows == 0 {
        return Err(BuildError::StructureParse(
            "manifest.json: contains no windows".to_string(),
        ));
    }
    Ok(ParsedManifest {
        probe,
        mode,
        n_windows,
    })
}

/// Render a chain-aware CHARMM `.psf` for an umbrella window (`PROA` for the
/// single amino-acid bead at index 0, `MEMB` for every lipid molecule).
pub fn render_umbrella_psf(particles: &[Particle], probe_name: &str) -> String {
    let lipid_ranges = molecule_ranges(&particles[1..]);
    let mut chains: Vec<(&str, std::ops::Range<usize>)> =
        Vec::with_capacity(1 + lipid_ranges.len());
    chains.push(("PROA", 0..1));
    for r in lipid_ranges {
        chains.push(("MEMB", (1 + r.start)..(1 + r.end)));
    }
    crate::psf::render_chains_ref(particles, &chains, probe_name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composition::Composition;
    use crate::membrane::{BuildOptions, Sizing};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

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
    fn parses_and_renders_lipid_protein_params_without_double_counting() {
        let ff = ForceField::isolf();
        let src = "\
# candidate trial from BO
[cg_ISOLF_parameters]
TRP 2.250 0.540 1.200

[cg_LJ_parameters]
TRP PHO 1.900 0.560
";
        let parsed = LipidProteinParams::parse(src).unwrap();
        let resolved = LipidProteinParams::resolve(
            &ff,
            &["TRP"],
            Some(&parsed),
            ProbeParamOverrides {
                auto_head_wca: true,
                ..ProbeParamOverrides::default()
            },
        )
        .unwrap();

        let itp = resolved.render_itp();
        assert!(itp.contains("[ cg_ISOLF_parameters ]"));
        assert!(itp.contains("   TRP    2.250    0.540    1.200"));
        assert!(itp.contains("[ cg_LJ_parameters ]"));
        assert!(itp.contains("   TRP   PHO    1.900    0.560    1.400"));
        // PHO is in LJ, so it must NOT also appear in WCA, whereas CHO/MID do.
        assert!(!resolved.wca.iter().any(|p| p.second == "PHO"));
        assert!(resolved.wca.iter().any(|p| p.second == "CHO"));
        assert!(resolved.wca.iter().any(|p| p.second == "MID"));
    }

    #[test]
    fn window_builder_preserves_lipid_count_and_centers_bilayer_at_zero() {
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
        assert!(mem.box_size[2] >= 20.0);

        let probe = Protein::single_bead("TRP", None).unwrap();
        let parts = build_window_particles(&mem, &probe, 1.5, &ff);
        assert_eq!(parts.len(), 1 + mem.particles.len());
        assert_eq!(parts[0].residue, "TRP");
        assert_eq!(parts[0].molecule_id, 1);
        assert_eq!(parts[0].atom_id, 1);
        assert!((parts[0].position[2] - 1.5).abs() < 1e-12);
        // Bilayer COM stays at z = 0.
        let mut msum = 0.0;
        let mut zsum = 0.0;
        for p in &parts[1..] {
            msum += p.mass;
            zsum += p.mass * p.position[2];
        }
        assert!((zsum / msum).abs() < 1e-10);
    }

    #[test]
    fn extracts_restr_cvs1_from_genesis_log() {
        let log = "\
INFO:       STEP            TIME       TOTAL_ENE   RESTRAINT_TOTAL   RESTRAINT001     RESTR_CVS001
------------ --------------- --------------- --------------- --------------- ---------------
INFO:          0          0.0000        -100.000             0.010             0.010           10.1414
INFO:        100          1.0000        -101.000             0.020             0.020            9.8500
INFO:        200          2.0000        -102.000             0.005             0.005           10.0707
";
        let rows = extract_cv_from_log(log, 100).unwrap();
        assert_eq!(rows, vec![(100, 9.85), (200, 10.0707)]);
        let dat = render_cv_dat(&rows);
        assert!(dat.contains("100       9.850000"));
        assert!(dat.contains("200      10.070700"));
    }

    #[test]
    fn demuxes_reus_logs_and_round_trips_manifest() {
        let log1 = "\
INFO: STEP RESTR_CVS001
INFO: 100  0.10
INFO: 200  5.10
";
        let log2 = "\
INFO: STEP RESTR_CVS001
INFO: 100  4.90
INFO: 200 -0.10
";
        let rem1 = "100 1\n200 2\n";
        let rem2 = "100 2\n200 1\n";
        let demux = demux_reus_cv(&[(log1, rem1), (log2, rem2)], 0).unwrap();
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
            mode: "reus",
            ensemble: "nvt",
            temperature_k: 323.15,
            k_kcal_mol_a2: 0.5,
            restraint_function: "DISTMASS",
            restraint_direction: "Z2",
            box_nm: [7.0, 7.0, 20.0],
            lipid_counts: &counts,
            z_windows_nm: &[0.0, 0.5, 1.0],
        };
        let json = manifest.to_json();
        assert!(json.contains("\"lipids\": {\"POPC\": 64, \"POPS\": 16}"));
        let parsed = parse_manifest(&json).unwrap();
        assert_eq!(parsed.probe, "TRP");
        assert_eq!(parsed.mode, "reus");
        assert_eq!(parsed.n_windows, 3);
    }
}
