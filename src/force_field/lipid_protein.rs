//! Tunable iSoLF lipid–protein non-bonded interaction parameters rendered into
//! a standalone `lipid_protein.itp` include file for GENESIS `RESIDCG`.
//!
//! GENESIS appends `[ cg_ISOLF_parameters ]`, `[ cg_WCA_parameters ]`, and
//! `[ cg_LJ_parameters ]` across `#include`d `.itp` files, so keeping
//! lipid–protein parameters in `lipid_protein.itp` leaves `isolf.itp` and
//! `aicg2p.itp` untouched and allows membrane+protein systems and umbrella
//! sampling / Bayesian-optimization campaigns to share the exact same
//! force-field layer.

use std::fmt::Write as _;

use crate::error::BuildError;
use crate::force_field::{ForceField, Interaction, LJ_CUTOFF_RATIO, POTENTIAL_MINIMUM_FACTOR};
use crate::protein::{Protein, params};

/// Default iSoLF well depth `epsilon` (kJ/mol) for an unparameterized amino-acid
/// bead before optimization.
pub const DEFAULT_PROBE_EPSILON: f64 = 1.5;

/// Default iSoLF polynomial-switch attraction width `omega` (nm).
pub const DEFAULT_PROBE_OMEGA: f64 = 1.0;

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

/// Inline CLI overrides applied to resolved `[ cg_ISOLF_parameters ]` entries.
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
    /// Duplicate names in `probe_aas` are deduplicated in first-seen order.
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
        let mut unique_aas: Vec<&str> = Vec::with_capacity(probe_aas.len());
        for &aa in probe_aas {
            if !unique_aas.contains(&aa) {
                unique_aas.push(aa);
            }
        }

        let mut out = base.cloned().unwrap_or_default();

        // Filter inherited entries from a master parameter file down to the
        // amino-acid types actually present in this system so unused residues
        // (e.g. an `ARG PHO` LJ pair in a `TRP`-only run) do not leak into
        // `lipid_protein.itp` or activate unused nonbonded kernels in GENESIS.
        if !unique_aas.is_empty() {
            let active = |name: &str| unique_aas.contains(&name);
            let active_pair = |a: &str, b: &str| {
                (active(a) && (active(b) || force_field.bead_type(b).is_some()))
                    || (active(b) && force_field.bead_type(a).is_some())
            };
            out.isolf.retain(|e| active(&e.name));
            out.wca.retain(|p| active_pair(&p.first, &p.second));
            out.lj.retain(|p| active_pair(&p.first, &p.second));
        }

        // Ensure every requested amino acid has an entry in `[ cg_ISOLF_parameters ]`
        // (unless the user explicitly omitted it while supplying WCA/LJ pairs for it).
        for &aa in &unique_aas {
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
            for &aa in &unique_aas {
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

    /// Resolve lipid–protein parameters for every unique amino-acid residue
    /// present across `proteins`.
    ///
    /// # Errors
    ///
    /// Propagates [`BuildError::InvalidSelection`] from [`Self::resolve`].
    pub fn for_proteins(
        force_field: &ForceField,
        proteins: &[&Protein],
        base: Option<&Self>,
        overrides: ProbeParamOverrides,
    ) -> Result<Self, BuildError> {
        let mut aas: Vec<&str> = Vec::new();
        for p in proteins {
            for r in p.unique_residues() {
                if !aas.contains(&r) {
                    aas.push(r);
                }
            }
        }
        Self::resolve(force_field, &aas, base, overrides)
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

#[cfg(test)]
mod tests {
    use super::*;

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
            &["TRP", "TRP"],
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
    fn resolves_parameters_across_multiple_proteins() {
        let ff = ForceField::isolf();
        let trp = Protein::single_bead("TRP", None).unwrap();
        let arg = Protein::single_bead("ARG", None).unwrap();
        let resolved = LipidProteinParams::for_proteins(
            &ff,
            &[&trp, &arg, &trp],
            None,
            ProbeParamOverrides {
                auto_head_wca: true,
                ..ProbeParamOverrides::default()
            },
        )
        .unwrap();
        assert_eq!(resolved.isolf.len(), 2);
        assert_eq!(resolved.isolf[0].name, "TRP");
        assert_eq!(resolved.isolf[1].name, "ARG");
        // 7 non-iSoLF head beads × 2 amino acids = 14 WCA pairs.
        assert_eq!(resolved.wca.len(), 14);
    }

    #[test]
    fn filters_master_table_to_active_probe_residues_only() {
        let ff = ForceField::isolf();
        let master = LipidProteinParams::parse(
            "\
[cg_ISOLF_parameters]
TRP 2.100 0.540 1.150
LEU 1.900 0.500 1.050
[cg_LJ_parameters]
ARG PHO 2.000 0.550
TRP PHO 1.800 0.560
",
        )
        .unwrap();
        let trp_only = LipidProteinParams::resolve(
            &ff,
            &["TRP"],
            Some(&master),
            ProbeParamOverrides {
                auto_head_wca: true,
                ..ProbeParamOverrides::default()
            },
        )
        .unwrap();
        assert_eq!(trp_only.isolf.len(), 1);
        assert_eq!(trp_only.isolf[0].name, "TRP");
        assert_eq!(trp_only.lj.len(), 1);
        assert_eq!(trp_only.lj[0].first, "TRP");
        assert!(trp_only.wca.iter().all(|p| p.first == "TRP"));
    }
}
