//! Parsers for CLI `KEY=VALUE` specifications, stage overrides, protein specs,
//! and umbrella probes.

use std::error::Error;
use std::fs;
use std::path::Path;

use isolf::composition::Composition;
use isolf::embed::{ProteinPlacement, ProteinSpec, ZPlacement};
use isolf::protein::Protein;
use isolf::protein::aicg2p::{
    AicgScale, ProteinModel, ProteinOptions, SafeDihedral, parse_charge_table, parse_selection,
};
use isolf::umbrella::{ProbeKind, UmbrellaProbe};
use isolf::vesicle::VesicleRadius;

/// Validated GENESIS `.inp` stage specs ready to write after the build finishes.
pub(crate) struct ControlSpecs {
    pub(crate) temperature: f64,
    pub(crate) ele_cutoff: Option<f64>,
    pub(crate) min: MinSpec,
    pub(crate) equil: DynSpec,
    pub(crate) prod: DynSpec,
}

/// Parsed minimization stage settings (before attaching the final box size).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct MinSpec {
    pub(crate) num_steps: usize,
    pub(crate) eneout_period: usize,
    pub(crate) crdout_period: usize,
    pub(crate) rstout_period: usize,
    pub(crate) nbupdate_period: usize,
}

impl Default for MinSpec {
    fn default() -> Self {
        Self {
            num_steps: 1000,
            eneout_period: 10,
            crdout_period: 1000,
            rstout_period: 1000,
            nbupdate_period: 10,
        }
    }
}

impl MinSpec {
    pub(crate) fn to_params(
        self,
        box_size: [f64; 3],
        temperature: f64,
        ele_cutoff: Option<f64>,
    ) -> isolf::inp::Minimization {
        let cg_ele_cutoff = ele_cutoff
            .map(isolf::inp::ele_cutoff_override)
            .or_else(|| isolf::inp::ele_cutoff_for_box(box_size));
        isolf::inp::Minimization {
            box_size,
            temperature,
            num_steps: self.num_steps,
            eneout_period: self.eneout_period,
            crdout_period: self.crdout_period,
            rstout_period: self.rstout_period,
            nbupdate_period: self.nbupdate_period,
            cg_ele_cutoff,
        }
    }
}

/// Parsed dynamics stage settings (`--equil` or `--prod`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct DynSpec {
    pub(crate) num_steps: usize,
    pub(crate) time_step: f64,
    pub(crate) temperature: Option<f64>,
    pub(crate) eneout_period: usize,
    pub(crate) crdout_period: usize,
    pub(crate) rstout_period: usize,
    pub(crate) stoptr_period: usize,
    pub(crate) nbupdate_period: usize,
    pub(crate) seed: Option<u64>,
}

impl DynSpec {
    pub(crate) fn default_equil() -> Self {
        Self {
            num_steps: 20_000,
            time_step: 0.01,
            temperature: None,
            eneout_period: 100,
            crdout_period: 100,
            rstout_period: 1000,
            stoptr_period: 10,
            nbupdate_period: 10,
            seed: None,
        }
    }

    pub(crate) fn default_prod() -> Self {
        Self {
            num_steps: 500_000,
            time_step: 0.02,
            temperature: None,
            eneout_period: 500,
            crdout_period: 500,
            rstout_period: 1000,
            stoptr_period: 10,
            nbupdate_period: 10,
            seed: None,
        }
    }

    pub(crate) fn to_params(
        self,
        default_temperature: f64,
        seed: u64,
        box_size: [f64; 3],
        ele_cutoff: Option<f64>,
    ) -> isolf::inp::Dynamics {
        let cg_ele_cutoff = ele_cutoff
            .map(isolf::inp::ele_cutoff_override)
            .or_else(|| isolf::inp::ele_cutoff_for_box(box_size));
        isolf::inp::Dynamics {
            temperature: self.temperature.unwrap_or(default_temperature),
            num_steps: self.num_steps,
            time_step: self.time_step,
            eneout_period: self.eneout_period,
            crdout_period: self.crdout_period,
            rstout_period: self.rstout_period,
            stoptr_period: self.stoptr_period,
            nbupdate_period: self.nbupdate_period,
            seed: self.seed.unwrap_or(seed),
            cg_ele_cutoff,
        }
    }
}

pub(crate) fn validate_padding(padding: f64) -> Result<(), Box<dyn Error>> {
    if !padding.is_finite() || padding < 0.0 {
        return Err(format!("--padding: '{padding}' must be a non-negative number of nm").into());
    }
    Ok(())
}

/// Resolve the two leaflet/shell compositions from symmetric `--lipids` or
/// paired (`--upper`/`--lower`, `--outer`/`--inner`) specs.
pub(crate) fn resolve_compositions(
    lipids: &[String],
    first: &[String],
    second: &[String],
    first_flag: &str,
    second_flag: &str,
) -> Result<(Composition, Composition), Box<dyn Error>> {
    let sym = join_specs(lipids);
    if !sym.is_empty() {
        let comp = Composition::parse(&sym)?;
        return Ok((comp.clone(), comp));
    }
    let a = join_specs(first);
    let b = join_specs(second);
    let a_opt = (!a.is_empty()).then_some(a.as_str());
    let b_opt = (!b.is_empty()).then_some(b.as_str());
    let first_spec = a_opt.or(b_opt).ok_or_else(|| {
        format!("missing lipid composition: pass --lipids <SPEC> (or {first_flag}/{second_flag})")
    })?;
    let second_spec = b_opt.or(a_opt).unwrap_or(first_spec);
    Ok((
        Composition::parse(first_spec)?,
        Composition::parse(second_spec)?,
    ))
}

/// Combine repeatable, comma-separated composition flags into one string.
pub(crate) fn join_specs(specs: &[String]) -> String {
    specs
        .iter()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join(",")
}

/// Parse a finite floating-point number (signed allowed), trimming whitespace.
pub(crate) fn parse_f64(value: &str, what: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("'{value}' is not a valid {what}").into())
}

/// Parse a non-negative finite floating-point number, trimming whitespace.
pub(crate) fn parse_nonneg_f64(value: &str, what: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v >= 0.0)
        .ok_or_else(|| format!("'{value}' is not a non-negative {what}").into())
}

/// Parse a single positive finite floating-point number, trimming whitespace.
pub(crate) fn parse_pos_f64(value: &str, unit: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v > 0.0)
        .ok_or_else(|| format!("'{value}' is not a positive {unit}").into())
}

/// Parse a positive integer (`>= 1`), trimming whitespace.
pub(crate) fn parse_pos_usize(value: &str, what: &str) -> Result<usize, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<usize>()
        .ok()
        .filter(|n| *n > 0)
        .ok_or_else(|| format!("'{value}' is not a positive {what}").into())
}

/// Parse a non-negative integer (`>= 0`), trimming whitespace.
pub(crate) fn parse_usize(value: &str, what: &str) -> Result<usize, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<usize>()
        .map_err(|_| format!("'{value}' is not a non-negative {what}").into())
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SpecTarget {
    /// `-p` on membrane/vesicle/tubule: accepts all placement and CG model keys.
    Embed,
    /// `--probe` on `umbrella setup` / `params`: accepts orientation (`tilt`,
    /// `rot`), `charge`, `name`, and CG model keys; rejects placement-only keys.
    Probe,
}

fn looks_like_structure_path(tok: &str) -> bool {
    Path::new(tok).exists() || tok.contains('/') || tok.contains('\\') || tok.contains('.')
}

/// Parse repeatable `-p, --protein <SPEC>[,KEY=VAL,...]` flags into [`ProteinSpec`]s.
pub(crate) fn parse_protein_specs(specs: &[String]) -> Result<Vec<ProteinSpec>, Box<dyn Error>> {
    let mut out = Vec::with_capacity(specs.len());
    for raw in specs {
        let (ps, _) = parse_single_protein_spec(SpecTarget::Embed, raw, None)?;
        out.push(ps);
    }
    Ok(out)
}

/// Parse `--probe <SPEC>` for `umbrella setup` and `params`:
/// - `all` or `aa:all` → all 20 canonical single-bead amino acids
/// - Bare or `aa:`-prefixed comma-separated AA list (`TRP`, `TRP,ARG`, `aa:TRP,ARG`)
///   when no `KEY=VAL` options are present and the first token is not a file path
/// - Single-bead AA with options (`TRP,charge=0`, `aa:HIS,charge=1.0,name=HIP`)
/// - PDB/mmCIF peptide or protein with optional orientation / CG model keys
///   (`walp23.pdb,tilt=90,rot=0,idr=1-3,charges=respac.dat`)
pub(crate) fn parse_umbrella_probes(
    spec: &str,
    charge_override: Option<f64>,
) -> Result<Vec<UmbrellaProbe>, Box<dyn Error>> {
    let trimmed = spec.trim();
    if trimmed.is_empty() {
        return Err("--probe: expected an amino acid (e.g. TRP, all) or PDB/mmCIF path".into());
    }

    let stripped_all = trimmed
        .strip_prefix("aa:")
        .or_else(|| trimmed.strip_prefix("AA:"))
        .unwrap_or(trimmed);
    if stripped_all.eq_ignore_ascii_case("all") {
        let mut out = Vec::with_capacity(20);
        for aa in Protein::all_amino_acid_names() {
            out.push(UmbrellaProbe {
                kind: ProbeKind::AminoAcid,
                protein: Protein::single_bead(aa, charge_override)?,
            });
        }
        return Ok(out);
    }

    // Check for a comma-separated list of bare amino acids (e.g. "TRP", "TRP,ARG",
    // "aa:TRP,ARG") without KEY=VAL options and not referencing a file path.
    let has_kv = trimmed.contains('=');
    let first_tok = trimmed.split(',').next().unwrap_or("").trim();

    if !has_kv && !looks_like_structure_path(first_tok) {
        let list_body = trimmed
            .strip_prefix("aa:")
            .or_else(|| trimmed.strip_prefix("AA:"))
            .unwrap_or(trimmed);
        let mut out: Vec<UmbrellaProbe> = Vec::new();
        for tok in list_body
            .split(',')
            .map(str::trim)
            .filter(|t| !t.is_empty())
        {
            let aa = tok
                .strip_prefix("aa:")
                .or_else(|| tok.strip_prefix("AA:"))
                .unwrap_or(tok);
            let p = Protein::single_bead(aa, charge_override)?;
            if !out.iter().any(|u| u.protein.name == p.name) {
                out.push(UmbrellaProbe {
                    kind: ProbeKind::AminoAcid,
                    protein: p,
                });
            }
        }
        if !out.is_empty() {
            return Ok(out);
        }
    }

    // Otherwise parse as a single spec (`TRP,charge=0`, `aa:TRP,charge=0`, or
    // `peptide.pdb,tilt=90,...`), rejecting placement-only keys that umbrella
    // setup controls automatically.
    let (mut ps, is_aa) = parse_single_protein_spec(SpecTarget::Probe, trimmed, charge_override)?;

    // Apply requested tilt (about x) and axial rotation (about z) around COM.
    let com = ps.protein.center_of_mass();
    if ps.placement.tilt_deg.abs() > 1e-12 {
        ps.protein
            .rotate_axis(com, [1.0, 0.0, 0.0], ps.placement.tilt_deg.to_radians());
    }
    if ps.placement.rot_deg.abs() > 1e-12 {
        ps.protein
            .rotate_axis(com, [0.0, 0.0, 1.0], ps.placement.rot_deg.to_radians());
    }

    Ok(vec![UmbrellaProbe {
        kind: if is_aa {
            ProbeKind::AminoAcid
        } else {
            ProbeKind::Protein
        },
        protein: ps.protein,
    }])
}

fn parse_single_protein_spec(
    target: SpecTarget,
    raw: &str,
    default_charge: Option<f64>,
) -> Result<(ProteinSpec, bool), Box<dyn Error>> {
    let flag = match target {
        SpecTarget::Embed => "--protein",
        SpecTarget::Probe => "--probe",
    };
    let mut tokens = raw.split(',').map(str::trim).filter(|t| !t.is_empty());
    let path_str = tokens
        .next()
        .ok_or_else(|| format!("{flag}: missing input structure or amino acid"))?;
    if path_str.contains('=') {
        return Err(format!(
            "{flag}: expected a PDB/mmCIF path or aa:<AA> before options, got '{path_str}'"
        )
        .into());
    }

    let mut placement = ProteinPlacement::default();
    let mut options = ProteinOptions::default();
    let mut name_override: Option<String> = None;
    let mut single_charge = default_charge;

    for tok in tokens {
        let (k, v) = tok
            .split_once('=')
            .ok_or_else(|| format!("{flag}: expected KEY=VAL option, got '{tok}'"))?;
        let key = k.trim().to_ascii_lowercase();
        let val = v.trim();
        if target == SpecTarget::Probe
            && matches!(
                key.as_str(),
                "z" | "zoff"
                    | "z_offset"
                    | "dz"
                    | "x"
                    | "y"
                    | "copies"
                    | "n"
                    | "clearance"
                    | "margin"
            )
        {
            return Err(format!(
                "--probe: key '{key}' is not used by umbrella probes \
                 (window z/x/y placement is automatic; valid keys: tilt, rot, charge, model, \
                 scale, safe-dihedral, idr, idr-local, idr-nonlocal, charges, contact-scale, name)"
            )
            .into());
        }
        match key.as_str() {
            "z" => {
                if val.eq_ignore_ascii_case("center") || val.eq_ignore_ascii_case("mid") {
                    placement.z_mode = ZPlacement::Center;
                } else if val.eq_ignore_ascii_case("keep") {
                    placement.z_mode = ZPlacement::Keep;
                } else {
                    placement.z_offset = parse_f64(val, "z offset in nm")?;
                }
            }
            "zoff" | "z_offset" | "dz" => {
                placement.z_offset = parse_f64(val, "z offset in nm")?;
            }
            "x" => placement.x = Some(parse_f64(val, "x coordinate in nm")?),
            "y" => placement.y = Some(parse_f64(val, "y coordinate in nm")?),
            "tilt" => placement.tilt_deg = parse_f64(val, "tilt angle in degrees")?,
            "rot" | "rotate" => placement.rot_deg = parse_f64(val, "rotation angle in degrees")?,
            "copies" | "n" => placement.copies = parse_pos_usize(val, "copy count")?,
            "clearance" | "margin" => {
                placement.clearance = parse_nonneg_f64(val, "clearance in nm")?;
            }
            "charge" | "q" => single_charge = Some(parse_f64(val, "charge in e")?),
            "model" | "ff" => {
                options.model = match val.to_ascii_lowercase().as_str() {
                    "aicg2+" | "aicg2p" | "aicg2" => ProteinModel::Aicg2Plus,
                    "clementi" | "ccgo" => ProteinModel::ClementiGo,
                    other => {
                        return Err(format!(
                            "{flag}: unknown model '{other}' (use aicg2+ or clementi)"
                        )
                        .into());
                    }
                };
            }
            "scale" => {
                options.scale = match val.to_ascii_lowercase().as_str() {
                    "general" | "1" => AicgScale::General,
                    "average" | "0" => AicgScale::Average,
                    other => {
                        return Err(format!(
                            "{flag}: unknown scale '{other}' (use general or average)"
                        )
                        .into());
                    }
                };
            }
            "safe-dihedral" | "safe_dihedral" | "safedih" => {
                options.safe_dihedral = match val.to_ascii_lowercase().as_str() {
                    "0" | "none" => SafeDihedral::None,
                    "1" | "cos2" => SafeDihedral::Cos2,
                    "3" | "sin3" => SafeDihedral::Sin3,
                    other => {
                        return Err(format!(
                            "{flag}: unknown safe-dihedral '{other}' (use 0, 1, or 3)"
                        )
                        .into());
                    }
                };
            }
            "idr" => {
                let norm = val.replace(['+', ';'], ",");
                let sel = parse_selection(&norm)?;
                options.idr_local = sel.clone();
                options.idr_nonlocal = sel;
            }
            "idr-local" | "idr_local" => {
                let norm = val.replace(['+', ';'], ",");
                options.idr_local = parse_selection(&norm)?;
            }
            "idr-nonlocal" | "idr_nonlocal" => {
                let norm = val.replace(['+', ';'], ",");
                options.idr_nonlocal = parse_selection(&norm)?;
            }
            "charges" | "respac" => {
                let content = fs::read_to_string(val)
                    .map_err(|e| format!("{flag}: failed to read charges file '{val}': {e}"))?;
                options.charges = parse_charge_table(&content)?;
            }
            "contact-scale" | "contact_scale" => {
                options.ccgo_contact_scale = parse_pos_f64(val, "contact scale")?;
            }
            "name" => name_override = Some(val.to_string()),
            other => {
                let valid = match target {
                    SpecTarget::Embed => {
                        "z, x, y, tilt, rot, copies, clearance, charge, model, scale, safe-dihedral, idr, idr-local, idr-nonlocal, charges, contact-scale, name"
                    }
                    SpecTarget::Probe => {
                        "tilt, rot, charge, model, scale, safe-dihedral, idr, idr-local, idr-nonlocal, charges, contact-scale, name"
                    }
                };
                return Err(format!("{flag}: unknown key '{other}' (use {valid})").into());
            }
        }
    }

    let explicit_aa = path_str
        .strip_prefix("aa:")
        .or_else(|| path_str.strip_prefix("AA:"));
    let bare_probe_aa = (target == SpecTarget::Probe
        && explicit_aa.is_none()
        && !looks_like_structure_path(path_str)
        && isolf::protein::params::canonical_aa(path_str).is_some_and(|c| c != "UNK"))
    .then_some(path_str);

    let (protein, is_aa) = if let Some(aa) = explicit_aa.or(bare_probe_aa) {
        let mut p = Protein::single_bead(aa, single_charge)?;
        if let Some(n) = name_override {
            p.name = isolf::protein::sanitize_mol_name(&n);
        }
        (p, true)
    } else {
        (
            Protein::from_file(Path::new(path_str), name_override.as_deref(), &options)?,
            false,
        )
    };
    Ok((ProteinSpec { protein, placement }, is_aa))
}

/// Walk non-empty comma-separated `KEY=VALUE` (or bare `VALUE`) tokens across
/// repeatable CLI flag occurrences. Returns `true` if any token was seen.
fn for_each_kv(
    specs: &[String],
    mut handle: impl FnMut(Option<&str>, &str) -> Result<(), Box<dyn Error>>,
) -> Result<bool, Box<dyn Error>> {
    let mut seen = false;
    for token in specs.iter().flat_map(|s| s.split(',')) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        seen = true;
        match token.split_once('=') {
            Some((k, v)) => handle(Some(k.trim()), v)?,
            None => handle(None, token)?,
        }
    }
    Ok(seen)
}

/// Parse `--count` occurrences into `(upper, lower)` lipid counts.
pub(crate) fn parse_count_spec(specs: &[String]) -> Result<Option<(usize, usize)>, Box<dyn Error>> {
    let (mut up, mut lo) = (None, None);
    let seen = for_each_kv(specs, |key, val| match key {
        None => {
            up = Some(parse_pos_usize(val, "lipid count")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("up") || k.eq_ignore_ascii_case("upper") => {
            up = Some(parse_pos_usize(val, "lipid count")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("lo") || k.eq_ignore_ascii_case("lower") => {
            lo = Some(parse_pos_usize(val, "lipid count")?);
            Ok(())
        }
        Some(other) => Err(format!("--count: unknown key '{other}' (use up or lo)").into()),
    })?;
    if !seen {
        return Ok(None);
    }
    match (up, lo) {
        (Some(u), Some(l)) => Ok(Some((u, l))),
        (Some(n), None) | (None, Some(n)) => Ok(Some((n, n))),
        (None, None) => Err("--count needs a value, e.g. `256` or `up=300,lo=200`".into()),
    }
}

/// Parse `--box` occurrences into `(x, y)` dimensions in nm.
pub(crate) fn parse_box_spec(specs: &[String]) -> Result<Option<(f64, f64)>, Box<dyn Error>> {
    let (mut x, mut y) = (None, None);
    let seen = for_each_kv(specs, |key, val| match key {
        None => {
            x = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("x") => {
            x = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("y") => {
            y = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(other) => Err(format!("--box: unknown key '{other}' (use x or y)").into()),
    })?;
    if !seen {
        return Ok(None);
    }
    match (x, y) {
        (Some(x), Some(y)) => Ok(Some((x, y))),
        (Some(s), None) | (None, Some(s)) => Ok(Some((s, s))),
        (None, None) => Err("--box needs a value, e.g. `20` or `x=30,y=15`".into()),
    }
}

/// Parse `--radius` occurrences into a [`VesicleRadius`].
pub(crate) fn parse_radius_spec(specs: &[String]) -> Result<Option<VesicleRadius>, Box<dyn Error>> {
    let (mut inner, mut outer) = (None, None);
    let seen = for_each_kv(specs, |key, val| match key {
        None => {
            outer = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("ri") || k.eq_ignore_ascii_case("inner") => {
            inner = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("ro") || k.eq_ignore_ascii_case("outer") => {
            outer = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(other) => Err(format!("--radius: unknown key '{other}' (use ro or ri)").into()),
    })?;
    if !seen {
        return Ok(None);
    }
    match (inner, outer) {
        (Some(inner), Some(outer)) => Ok(Some(VesicleRadius::Both { inner, outer })),
        (Some(inner), None) => Ok(Some(VesicleRadius::Inner(inner))),
        (None, Some(outer)) => Ok(Some(VesicleRadius::Outer(outer))),
        (None, None) => Err("--radius needs a value, e.g. `20` or `ri=14,ro=20`".into()),
    }
}

/// Parse `--min` override tokens into a [`MinSpec`].
pub(crate) fn parse_min_spec(specs: &[String]) -> Result<MinSpec, Box<dyn Error>> {
    let mut out = MinSpec::default();
    for_each_kv(specs, |key, val| match key {
        None => {
            out.num_steps = parse_pos_usize(val, "step count")?;
            Ok(())
        }
        Some(k) => match k.to_ascii_lowercase().as_str() {
            "steps" | "nsteps" | "num_steps" => {
                out.num_steps = parse_pos_usize(val, "step count")?;
                Ok(())
            }
            "eneout" | "eneout_period" => {
                out.eneout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "crdout" | "crdout_period" => {
                out.crdout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "rstout" | "rstout_period" => {
                out.rstout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "nbupdate" | "nbupdate_period" => {
                out.nbupdate_period = parse_usize(val, "period")?;
                Ok(())
            }
            other => Err(format!(
                "--min: unknown key '{other}' (use steps, eneout, crdout, rstout, nbupdate)"
            )
            .into()),
        },
    })?;
    Ok(out)
}

/// Parse `--equil` or `--prod` override tokens into a [`DynSpec`].
pub(crate) fn parse_dyn_spec(
    flag: &str,
    specs: &[String],
    mut out: DynSpec,
) -> Result<DynSpec, Box<dyn Error>> {
    for_each_kv(specs, |key, val| {
        match key {
        None => {
            out.num_steps = parse_pos_usize(val, "step count")?;
            Ok(())
        }
        Some(k) => match k.to_ascii_lowercase().as_str() {
            "steps" | "nsteps" | "num_steps" => {
                out.num_steps = parse_pos_usize(val, "step count")?;
                Ok(())
            }
            "dt" | "timestep" | "time_step" => {
                out.time_step = parse_pos_f64(val, "time step")?;
                Ok(())
            }
            "temp" | "temperature" => {
                out.temperature = Some(parse_pos_f64(val, "temperature in K")?);
                Ok(())
            }
            "eneout" | "eneout_period" => {
                out.eneout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "crdout" | "crdout_period" => {
                out.crdout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "rstout" | "rstout_period" => {
                out.rstout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "stoptr" | "stoptr_period" => {
                out.stoptr_period = parse_usize(val, "period")?;
                Ok(())
            }
            "nbupdate" | "nbupdate_period" => {
                out.nbupdate_period = parse_usize(val, "period")?;
                Ok(())
            }
            "seed" | "iseed" => {
                let s = val
                    .trim()
                    .parse::<u64>()
                    .map_err(|_| format!("'{val}' is not a valid seed"))?;
                out.seed = Some(s);
                Ok(())
            }
            other => Err(format!(
                "{flag}: unknown key '{other}' (use steps, dt, temp, eneout, crdout, rstout, stoptr, nbupdate, seed)"
            )
            .into()),
        },
    }
    })?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_geometry_specs() {
        assert_eq!(parse_count_spec(&["256".into()]).unwrap(), Some((256, 256)));
        assert_eq!(
            parse_count_spec(&["up=300,lo=200".into()]).unwrap(),
            Some((300, 200))
        );
        assert!(parse_count_spec(&["foo=10".into()]).is_err());

        assert_eq!(parse_box_spec(&["15".into()]).unwrap(), Some((15.0, 15.0)));
        assert_eq!(
            parse_box_spec(&["x=10".into(), "y=20".into()]).unwrap(),
            Some((10.0, 20.0))
        );
        assert!(parse_box_spec(&["z=10".into()]).is_err());

        assert_eq!(
            parse_radius_spec(&["20".into()]).unwrap(),
            Some(VesicleRadius::Outer(20.0))
        );
        assert_eq!(
            parse_radius_spec(&["ri=14,ro=20".into()]).unwrap(),
            Some(VesicleRadius::Both {
                inner: 14.0,
                outer: 20.0
            })
        );
        assert!(parse_radius_spec(&["r=20".into()]).is_err());
    }

    #[test]
    fn parses_stage_override_specs() {
        let m = parse_min_spec(&["2000".into(), "eneout=20,crdout=500".into()]).unwrap();
        assert_eq!(m.num_steps, 2000);
        assert_eq!(m.eneout_period, 20);
        assert_eq!(m.crdout_period, 500);
        assert_eq!(m.nbupdate_period, 10);
        assert!(parse_min_spec(&["dt=0.01".into()]).is_err());

        let d = parse_dyn_spec(
            "--equil",
            &["steps=50000,dt=0.015,temp=310,seed=7".into()],
            DynSpec::default_equil(),
        )
        .unwrap();
        assert_eq!(d.num_steps, 50_000);
        assert_eq!(d.time_step, 0.015);
        assert_eq!(d.temperature, Some(310.0));
        assert_eq!(d.seed, Some(7));
        assert!(parse_dyn_spec("--prod", &["bad=1".into()], DynSpec::default_prod()).is_err());
    }

    #[test]
    fn parses_umbrella_probe_variants() {
        let one = parse_umbrella_probes("TRP", None).unwrap();
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].kind, ProbeKind::AminoAcid);
        assert_eq!(one[0].protein.name, "TRP");

        let two = parse_umbrella_probes("aa:TRP,ARG", None).unwrap();
        assert_eq!(two.len(), 2);
        assert_eq!(two[1].protein.name, "ARG");

        let all = parse_umbrella_probes("all", None).unwrap();
        assert_eq!(all.len(), 20);

        let custom = parse_umbrella_probes("aa:HIS,charge=1.0,name=HIP", None).unwrap();
        assert_eq!(custom[0].kind, ProbeKind::AminoAcid);
        assert_eq!(custom[0].protein.name, "HIP");
        assert!((custom[0].protein.beads[0].charge - 1.0).abs() < 1e-12);

        // Bare AA + KEY=VAL works without requiring the `aa:` prefix.
        let bare_kv = parse_umbrella_probes("TRP,charge=0.0", None).unwrap();
        assert_eq!(bare_kv[0].kind, ProbeKind::AminoAcid);
        assert_eq!(bare_kv[0].protein.name, "TRP");
        assert!((bare_kv[0].protein.beads[0].charge - 0.0).abs() < 1e-12);

        // Placement-only keys are rejected instead of being silently ignored.
        assert!(parse_umbrella_probes("TRP,z=1.5", None).is_err());
        assert!(parse_umbrella_probes("peptide.pdb,copies=2", None).is_err());
        assert!(parse_umbrella_probes("peptide.pdb,clearance=0.3", None).is_err());
    }
}
