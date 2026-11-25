//! GENESIS RESIDCG control files (`.inp`) for a three-stage iSoLF run: energy
//! minimization, then dynamics (NVT or NPT, chosen by geometry) for an
//! equilibration stage and a production stage.
//!
//! The minimization file carries an explicit `box_size`: GENESIS requires it for
//! the first (restart-free) run and does not read the box from the `.gro`. Later
//! stages read the box from the preceding restart. The `.gro` box is in nm and
//! GENESIS expects Ångström, so [`minimization`] converts it.
//!
//! `nsteps` must be an exact multiple of every output/update period, otherwise
//! the run would stop mid-interval; the generators enforce this and return an
//! [`InpError`] on a mismatch.

use std::fmt::{self, Write as _};

/// nm → Ångström, the unit GENESIS uses internally and in `box_size`.
const ANGSTROM_PER_NM: f64 = 10.0;

/// GENESIS `RESIDCG` default Debye–Hückel electrostatic cutoff (`cg_cutoffdist_ele`, Å).
const DEFAULT_CG_CUTOFF_ELE_A: f64 = 52.0;

/// GENESIS `RESIDCG` default Verlet skin (`cg_pairlistdist_ele - cg_cutoffdist_ele`, Å).
const DEFAULT_CG_ELE_SKIN_A: f64 = 5.0;

/// Minimum box dimension (Å) above which GENESIS's built-in 52 / 57 Å
/// electrostatic cutoffs fit comfortably with room for NPT box compression.
const SAFE_DEFAULT_ELE_BOX_A: f64 = 120.0;

/// Compute `(cg_cutoffdist_ele, cg_pairlistdist_ele)` in Ångström when the
/// simulation box is too small for GENESIS `RESIDCG`'s built-in 52 / 57 Å
/// electrostatic defaults (`box > 2 * cutoff`).
///
/// Returns `None` when every box axis is at least `12.0 nm` (120 Å), letting
/// GENESIS keep its built-in defaults untouched.
pub fn ele_cutoff_for_box(box_nm: [f64; 3]) -> Option<(f64, f64)> {
    let min_box_a = box_nm[0].min(box_nm[1]).min(box_nm[2]) * ANGSTROM_PER_NM;
    if min_box_a >= SAFE_DEFAULT_ELE_BOX_A {
        return None;
    }
    let cutoff = (0.42 * min_box_a)
        .floor()
        .clamp(1.0, DEFAULT_CG_CUTOFF_ELE_A);
    let max_pairlist = (0.48 * min_box_a).floor().max(cutoff + 1.0);
    let pairlist = (cutoff + DEFAULT_CG_ELE_SKIN_A).min(max_pairlist);
    Some((cutoff, pairlist))
}

/// Pair an explicit user-supplied `cg_cutoffdist_ele` (Å) with its matching
/// Verlet pairlist distance (`cutoff + 5.0 Å`).
pub fn ele_cutoff_override(cutoff_angstrom: f64) -> (f64, f64) {
    (cutoff_angstrom, cutoff_angstrom + DEFAULT_CG_ELE_SKIN_A)
}

/// Error from generating a control file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InpError {
    /// `nsteps` is not an exact multiple of the named output/update period.
    NotMultipleOf(&'static str),
}

impl fmt::Display for InpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotMultipleOf(name) => write!(f, "nsteps must be a multiple of {name}"),
        }
    }
}

impl std::error::Error for InpError {}

/// Require `num_steps` to be an exact multiple of `period` (a period of 0
/// disables that output, so it carries no constraint).
fn require_multiple(num_steps: usize, period: usize, name: &'static str) -> Result<(), InpError> {
    if period != 0 && !num_steps.is_multiple_of(period) {
        return Err(InpError::NotMultipleOf(name));
    }
    Ok(())
}

/// Parameters for the energy-minimization (steepest descent) control file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Minimization {
    /// Periodic box (nm); written to `box_size_*` in Ångström.
    pub box_size: [f64; 3],
    /// Implicit-solvent temperature (K), `cg_sol_temperature`.
    pub temperature: f64,
    /// Number of minimization steps.
    pub num_steps: usize,
    /// Energy output period.
    pub eneout_period: usize,
    /// Coordinate output period.
    pub crdout_period: usize,
    /// Restart output period.
    pub rstout_period: usize,
    /// Neighbour-list update period.
    pub nbupdate_period: usize,
    /// Optional `(cg_cutoffdist_ele, cg_pairlistdist_ele)` override in Ångström
    /// for small boxes where GENESIS's default 52 / 57 Å exceeds `box / 2`.
    pub cg_ele_cutoff: Option<(f64, f64)>,
}

/// Parameters for an NVT or NPT dynamics control file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dynamics {
    /// Temperature (K); written to both `temperature` and `cg_sol_temperature`.
    pub temperature: f64,
    /// Number of dynamics steps.
    pub num_steps: usize,
    /// Integration time step.
    pub time_step: f64,
    /// Energy output period.
    pub eneout_period: usize,
    /// Coordinate output period.
    pub crdout_period: usize,
    /// Restart output period.
    pub rstout_period: usize,
    /// Translational-motion removal period.
    pub stoptr_period: usize,
    /// Neighbour-list update period.
    pub nbupdate_period: usize,
    /// Langevin thermostat seed (`iseed`).
    pub seed: u64,
    /// Optional `(cg_cutoffdist_ele, cg_pairlistdist_ele)` override in Ångström
    /// for small boxes where GENESIS's default 52 / 57 Å exceeds `box / 2`.
    pub cg_ele_cutoff: Option<(f64, f64)>,
}

/// Single-window signed z-distance umbrella restraint (`DISTMASS`/`DIST` with
/// `direction1 = Z2`: `U = k * (dz - z0)^2`, `dz = z1 - z2` in Å).
#[derive(Debug, Clone, PartialEq)]
pub struct ZRestraint<'a> {
    /// GENESIS atom selection for group 1 (the probe, e.g. `"ano:1"` or `"ano:1-23"`).
    pub group1: &'a str,
    /// GENESIS atom selection for group 2 (the bilayer, e.g. `"rnam:POPC"`).
    pub group2: &'a str,
    /// Restraint function (`"DISTMASS"` for COM or `"DIST"` for geometric center).
    pub function: &'a str,
    /// Direction (`"Z2"` for signed `z1 - z2` along the membrane normal).
    pub direction: &'a str,
    /// Harmonic force constant `k` in kcal·mol⁻¹·Å⁻² (`U = k * (dz - z0)^2`).
    pub constant: f64,
    /// Target signed z-distance `z0 = z1 - z2` in Ångström.
    pub reference_angstrom: f64,
}

/// Multi-window Replica-Exchange Umbrella Sampling (`[REMD]` `type1 = RESTRAINT`)
/// settings for GENESIS `atdyn`.
#[derive(Debug, Clone, PartialEq)]
pub struct ReusConfig<'a> {
    /// Replica-exchange attempt interval in MD steps (must divide `nsteps`).
    pub exchange_period: usize,
    /// RNG seed for replica-exchange Monte Carlo (`[REMD] iseed`).
    pub iseed: u64,
    /// GENESIS atom selection for group 1 (the probe).
    pub group1: &'a str,
    /// GENESIS atom selection for group 2 (the bilayer).
    pub group2: &'a str,
    /// Restraint function (`"DISTMASS"` or `"DIST"`).
    pub function: &'a str,
    /// Direction (`"Z2"`).
    pub direction: &'a str,
    /// Per-replica force constants `k` in kcal·mol⁻¹·Å⁻² (`U = k * (dz - z0)^2`).
    pub constants: &'a [f64],
    /// Per-replica target signed z-distances `z0` in Ångström.
    pub references_angstrom: &'a [f64],
}

/// Configuration for a GENESIS `mbar_analysis` 1D umbrella-sampling PMF file.
#[derive(Debug, Clone, PartialEq)]
pub struct MbarConfig<'a> {
    /// CV trajectory file pattern with `{}` replica placeholder (e.g. `"./cv/{}.dat"`).
    pub cvfile: &'a str,
    /// Output relative free-energy file path.
    pub fenefile: &'a str,
    /// Output PMF file path.
    pub pmffile: &'a str,
    /// Optional per-window weight file pattern (with `{}`). GENESIS
    /// `mbar_analysis` only supports writing weights when `nblocks == 1`.
    pub weightfile: Option<&'a str>,
    /// Simulation and target temperature in K.
    pub temperature: f64,
    /// Number of blocks for block-averaging error estimation.
    pub nblocks: usize,
    /// Lower bound of the PMF grid in Å.
    pub grid_min_angstrom: f64,
    /// Upper bound of the PMF grid in Å.
    pub grid_max_angstrom: f64,
    /// Number of grid points (bins + 1).
    pub grid_points: usize,
    /// Per-window force constants `k` in kcal·mol⁻¹·Å⁻² (matches GENESIS `U = k*(dz - z0)^2`).
    pub constants: &'a [f64],
    /// Per-window target signed z-distances `z0` in Ångström.
    pub references_angstrom: &'a [f64],
}

/// Energy-minimization control file reading the `<structure>.{top,gro}` pair.
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output or update period.
pub fn minimization(structure: &str, params: &Minimization) -> Result<String, InpError> {
    minimization_with_restraint(
        &format!("./{structure}.top"),
        &format!("./{structure}.gro"),
        "min",
        params,
        None,
    )
}

/// Energy-minimization control file with custom topology/coordinate paths,
/// output basename, and an optional umbrella z-restraint.
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output or update period.
pub fn minimization_with_restraint(
    top_path: &str,
    gro_path: &str,
    basename: &str,
    params: &Minimization,
    restraint: Option<&ZRestraint<'_>>,
) -> Result<String, InpError> {
    let Minimization {
        box_size,
        temperature,
        num_steps,
        eneout_period,
        crdout_period,
        rstout_period,
        nbupdate_period,
        cg_ele_cutoff,
    } = *params;
    require_multiple(num_steps, eneout_period, "eneout_period")?;
    require_multiple(num_steps, crdout_period, "crdout_period")?;
    require_multiple(num_steps, rstout_period, "rstout_period")?;
    require_multiple(num_steps, nbupdate_period, "nbupdate_period")?;

    let [x, y, z] = box_size.map(|nm| nm * ANGSTROM_PER_NM);
    let ele_cutoff_lines = format_ele_cutoff(cg_ele_cutoff);
    let mut s = format!(
        "\
[INPUT]
grotopfile            = {top_path}
grocrdfile            = {gro_path}

[OUTPUT]
rstfile               = ./{basename}.rst
dcdfile               = ./{basename}.dcd
pdbfile               = ./{basename}.pdb

[ENERGY]
forcefield            = RESIDCG
electrostatic         = CUTOFF
cg_sol_temperature    = {temperature}
cg_sol_ionic_strength = 0.15
{ele_cutoff_lines}
[MINIMIZE]
method                = SD
nsteps                = {num_steps}
eneout_period         = {eneout_period}
crdout_period         = {crdout_period}
rstout_period         = {rstout_period}
nbupdate_period       = {nbupdate_period}
check_structure       = NO

[BOUNDARY]
type                  = PBC
box_size_x            = {x:.4}
box_size_y            = {y:.4}
box_size_z            = {z:.4}
"
    );
    if let Some(r) = restraint {
        append_single_restraint(&mut s, r);
    }
    Ok(s)
}

/// Thermostat ensemble for a dynamics stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ensemble {
    /// Constant volume (vesicles, tubules, and fixed-area umbrella sampling).
    Nvt,
    /// Constant pressure, semi-isotropic and Z-fixed (flat membranes).
    Npt,
}

/// An NVT control file over the `<structure>.{top,gro}` pair, reading `restart_in`
/// and writing `<basename>.{rst,dcd,pdb}`.
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output or update period.
pub fn nvt(
    structure: &str,
    restart_in: &str,
    basename: &str,
    params: &Dynamics,
) -> Result<String, InpError> {
    dynamics_with_restraint(
        Ensemble::Nvt,
        &format!("./{structure}.top"),
        &format!("./{structure}.gro"),
        &format!("./{restart_in}"),
        basename,
        params,
        None,
    )
}

/// An NPT (semi-isotropic, Z-fixed) control file over the `<structure>.{top,gro}`
/// pair, reading `restart_in` and writing `<basename>.{rst,dcd,pdb}`. The barostat
/// relaxes the leaflet to its tensionless area per lipid.
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output or update period.
pub fn npt(
    structure: &str,
    restart_in: &str,
    basename: &str,
    params: &Dynamics,
) -> Result<String, InpError> {
    dynamics_with_restraint(
        Ensemble::Npt,
        &format!("./{structure}.top"),
        &format!("./{structure}.gro"),
        &format!("./{restart_in}"),
        basename,
        params,
        None,
    )
}

/// Render a dynamics control file with explicit file paths and an optional
/// umbrella z-restraint (`[SELECTION]` + `[RESTRAINTS]`).
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output or update period.
pub fn dynamics_with_restraint(
    ensemble: Ensemble,
    top_path: &str,
    gro_path: &str,
    restart_in_path: &str,
    basename: &str,
    params: &Dynamics,
    restraint: Option<&ZRestraint<'_>>,
) -> Result<String, InpError> {
    validate_dynamics(params)?;
    let body = dynamics_body(ensemble, params);
    let mut s = format!(
        "\
[INPUT]
grotopfile            = {top_path}
grocrdfile            = {gro_path}
rstfile               = {restart_in_path}

[OUTPUT]
rstfile               = ./{basename}.rst
dcdfile               = ./{basename}.dcd
pdbfile               = ./{basename}.pdb

{body}"
    );
    if let Some(r) = restraint {
        append_single_restraint(&mut s, r);
    }
    Ok(s)
}

/// Render a GENESIS Replica-Exchange Umbrella Sampling (`[REMD]` `type1 = RESTRAINT`)
/// control file for `atdyn`.
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output/update/exchange period.
pub fn reus(
    ensemble: Ensemble,
    top_path: &str,
    gro_path: &str,
    restart_in_pattern: &str,
    basename: &str,
    params: &Dynamics,
    reus_cfg: &ReusConfig<'_>,
) -> Result<String, InpError> {
    validate_dynamics(params)?;
    require_multiple(
        params.num_steps,
        reus_cfg.exchange_period,
        "exchange_period",
    )?;
    let nreplica = reus_cfg.references_angstrom.len();
    let body = dynamics_body(ensemble, params);
    let mut s = format!(
        "\
[INPUT]
grotopfile            = {top_path}
grocrdfile            = {gro_path}
rstfile               = {restart_in_pattern}

[OUTPUT]
logfile               = ./{basename}_{{}}.log
rstfile               = ./{basename}_{{}}.rst
dcdfile               = ./{basename}_{{}}.dcd
remfile               = ./{basename}_{{}}.rem

[REMD]
dimension             = 1
exchange_period       = {exchange_period}
iseed                 = {remd_seed}
type1                 = RESTRAINT
nreplica1             = {nreplica}
cyclic_params1        = NO
rest_function1        = 1

{body}
[SELECTION]
group1                = {group1}
group2                = {group2}

[RESTRAINTS]
nfunctions            = 1
function1             = {function}
direction1            = {direction}
select_index1         = 1 2
",
        exchange_period = reus_cfg.exchange_period,
        remd_seed = reus_cfg.iseed,
        group1 = reus_cfg.group1,
        group2 = reus_cfg.group2,
        function = reus_cfg.function,
        direction = reus_cfg.direction,
    );
    write_wrapped_floats(&mut s, "constant1             =", reus_cfg.constants);
    write_wrapped_floats(
        &mut s,
        "reference1            =",
        reus_cfg.references_angstrom,
    );
    Ok(s)
}

/// Shared `[ENERGY]`, `[DYNAMICS]`, `[CONSTRAINTS]`, `[ENSEMBLE]`, and
/// `[BOUNDARY]` sections for standard dynamics and REUS control files.
fn dynamics_body(ensemble: Ensemble, params: &Dynamics) -> String {
    let Dynamics {
        temperature,
        num_steps,
        time_step,
        eneout_period,
        crdout_period,
        rstout_period,
        stoptr_period,
        nbupdate_period,
        seed,
        cg_ele_cutoff,
    } = *params;
    let ensemble = ensemble_block(ensemble, temperature);
    let ele_cutoff_lines = format_ele_cutoff(cg_ele_cutoff);
    format!(
        "\
[ENERGY]
forcefield            = RESIDCG
electrostatic         = CUTOFF
cg_sol_temperature    = {temperature}
cg_sol_ionic_strength = 0.15
{ele_cutoff_lines}
[DYNAMICS]
integrator            = VVER_CG
timestep              = {time_step}
nsteps                = {num_steps}
eneout_period         = {eneout_period}
crdout_period         = {crdout_period}
rstout_period         = {rstout_period}
stoptr_period         = {stoptr_period}
nbupdate_period       = {nbupdate_period}
iseed                 = {seed}

[CONSTRAINTS]
rigid_bond            = NO

{ensemble}
[BOUNDARY]
type                  = PBC
"
    )
}

/// Render a GENESIS `mbar_analysis` control file for 1D umbrella sampling
/// along the bilayer normal (`U = k * (dz - z0)^2`, `z` in Å, output in kcal/mol).
pub fn mbar(cfg: &MbarConfig<'_>) -> String {
    let nreplica = cfg.references_angstrom.len();
    let weight_line = match cfg.weightfile {
        Some(w) if cfg.nblocks == 1 => format!("weightfile         = {w}\n"),
        _ => String::new(),
    };
    let mut s = format!(
        "\
[INPUT]
cvfile             = {cvfile}

[OUTPUT]
fenefile           = {fenefile}
pmffile            = {pmffile}
{weight_line}
[MBAR]
nreplica           = {nreplica}
input_type         = US
dimension          = 1
nblocks            = {nblocks}
self_iteration     = 5
newton_iteration   = 100
temperature        = {temp}
target_temperature = {temp}
tolerance          = 1.0E-08
rest_function1     = 1
grids1             = {gmin:.4} {gmax:.4} {gpts}
output_unit        = kcal/mol

[RESTRAINTS]
",
        cvfile = cfg.cvfile,
        fenefile = cfg.fenefile,
        pmffile = cfg.pmffile,
        nblocks = cfg.nblocks,
        temp = cfg.temperature,
        gmin = cfg.grid_min_angstrom,
        gmax = cfg.grid_max_angstrom,
        gpts = cfg.grid_points,
    );
    write_wrapped_floats(&mut s, "constant1          =", cfg.constants);
    write_wrapped_floats(&mut s, "reference1         =", cfg.references_angstrom);
    s.push_str("is_periodic1       = NO\n");
    s
}

fn format_ele_cutoff(cg_ele_cutoff: Option<(f64, f64)>) -> String {
    match cg_ele_cutoff {
        Some((rc, rp)) => {
            format!("cg_cutoffdist_ele     = {rc:.1}\ncg_pairlistdist_ele   = {rp:.1}\n")
        }
        None => String::new(),
    }
}

fn append_single_restraint(s: &mut String, r: &ZRestraint<'_>) {
    let _ = write!(
        s,
        "\n\
[SELECTION]
group1                = {group1}
group2                = {group2}

[RESTRAINTS]
nfunctions            = 1
function1             = {function}
direction1            = {direction}
select_index1         = 1 2
constant1             = {constant:.4}
reference1            = {reference:.4}
",
        group1 = r.group1,
        group2 = r.group2,
        function = r.function,
        direction = r.direction,
        constant = r.constant,
        reference = r.reference_angstrom,
    );
}

/// Write a GENESIS key followed by a space-separated float array, wrapping long
/// arrays across continuation lines (`\`) at 10 values per row.
fn write_wrapped_floats(s: &mut String, key_prefix: &str, values: &[f64]) {
    if values.len() <= 10 {
        s.push_str(key_prefix);
        for &v in values {
            let _ = write!(s, " {v:.4}");
        }
        s.push('\n');
        return;
    }
    let _ = writeln!(s, "{key_prefix} \\");
    for (chunk_idx, chunk) in values.chunks(10).enumerate() {
        let is_last = (chunk_idx + 1) * 10 >= values.len();
        s.push(' ');
        for &v in chunk {
            let _ = write!(s, " {v:>8.4}");
        }
        if !is_last {
            s.push_str(" \\");
        }
        s.push('\n');
    }
}

/// The `[ENSEMBLE]` block for the chosen ensemble at `temperature`.
fn ensemble_block(ensemble: Ensemble, temperature: f64) -> String {
    match ensemble {
        Ensemble::Nvt => format!(
            "\
[ENSEMBLE]
ensemble              = NVT
tpcontrol             = LANGEVIN
temperature           = {temperature}
pressure              = 0.00
gamma_t               = 0.01
"
        ),
        Ensemble::Npt => format!(
            "\
[ENSEMBLE]
ensemble              = NPT
tpcontrol             = LANGEVIN
temperature           = {temperature}
pressure              = 0.00
gamma_t               = 0.01
gamma_p               = 0.01
isotropy              = Z-FIXED-SEMI-ISO
"
        ),
    }
}

/// Check that an NVT/NPT step count is a multiple of every period it writes.
fn validate_dynamics(params: &Dynamics) -> Result<(), InpError> {
    require_multiple(params.num_steps, params.eneout_period, "eneout_period")?;
    require_multiple(params.num_steps, params.crdout_period, "crdout_period")?;
    require_multiple(params.num_steps, params.rstout_period, "rstout_period")?;
    require_multiple(params.num_steps, params.stoptr_period, "stoptr_period")?;
    require_multiple(params.num_steps, params.nbupdate_period, "nbupdate_period")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dynamics() -> Dynamics {
        Dynamics {
            temperature: 323.15,
            num_steps: 1000,
            time_step: 0.04,
            eneout_period: 250,
            crdout_period: 250,
            rstout_period: 1000,
            stoptr_period: 10,
            nbupdate_period: 10,
            seed: 42,
            cg_ele_cutoff: None,
        }
    }

    fn minimization_params() -> Minimization {
        Minimization {
            box_size: [2.5, 3.0, 5.0],
            temperature: 323.15,
            num_steps: 1000,
            eneout_period: 10,
            crdout_period: 1000,
            rstout_period: 1000,
            nbupdate_period: 10,
            cg_ele_cutoff: None,
        }
    }

    #[test]
    fn minimization_writes_box_in_angstrom_and_parameters() {
        let inp = minimization("membrane", &minimization_params()).unwrap();
        assert!(inp.contains("grotopfile            = ./membrane.top"));
        assert!(inp.contains("grocrdfile            = ./membrane.gro"));
        // nm × 10 → Å.
        assert!(inp.contains("box_size_x            = 25.0000"));
        assert!(inp.contains("box_size_y            = 30.0000"));
        assert!(inp.contains("box_size_z            = 50.0000"));
        assert!(inp.contains("cg_sol_temperature    = 323.15"));
        assert!(inp.contains("nsteps                = 1000"));
        assert!(inp.contains("eneout_period         = 10"));
        assert!(inp.contains("nbupdate_period       = 10"));
        assert!(inp.contains("check_structure       = NO"));
        assert!(inp.contains("method                = SD"));
        // Minimization has no time step or motion-removal period.
        assert!(!inp.contains("timestep"));
        assert!(!inp.contains("stoptr_period"));
    }

    #[test]
    fn dynamics_writes_periods_temperature_seed_and_restart_chain() {
        let nvt_inp = nvt(
            "membrane",
            "min.rst",
            "nvt",
            &Dynamics {
                temperature: 310.0,
                seed: 111,
                num_steps: 1000,
                ..dynamics()
            },
        )
        .unwrap();
        assert!(nvt_inp.contains("grotopfile            = ./membrane.top"));
        assert!(nvt_inp.contains("ensemble              = NVT"));
        assert!(nvt_inp.contains("cg_sol_temperature    = 310"));
        assert!(nvt_inp.contains("temperature           = 310"));
        assert!(nvt_inp.contains("iseed                 = 111"));
        assert!(nvt_inp.contains("nsteps                = 1000"));
        assert!(nvt_inp.contains("stoptr_period         = 10"));
        assert!(nvt_inp.contains("nbupdate_period       = 10"));
        assert!(nvt_inp.contains("rstfile               = ./min.rst")); // reads the given restart
        assert!(nvt_inp.contains("rstfile               = ./nvt.rst")); // writes the given basename

        let npt_inp = npt(
            "membrane",
            "min.rst",
            "npt",
            &Dynamics {
                seed: 222,
                num_steps: 25000,
                ..dynamics()
            },
        )
        .unwrap();
        assert!(npt_inp.contains("ensemble              = NPT"));
        assert!(npt_inp.contains("isotropy              = Z-FIXED-SEMI-ISO"));
        assert!(npt_inp.contains("iseed                 = 222"));
        assert!(npt_inp.contains("nsteps                = 25000"));
        assert!(npt_inp.contains("stoptr_period         = 10"));

        // A production stage is the same generator chained off the equilibration;
        // the structure name flows into grotopfile/grocrdfile.
        let pro = npt("vesicle", "npt.rst", "pro", &dynamics()).unwrap();
        assert!(pro.contains("grotopfile            = ./vesicle.top"));
        assert!(pro.contains("rstfile               = ./npt.rst")); // reads npt restart
        assert!(pro.contains("rstfile               = ./pro.rst")); // writes pro
        assert!(pro.contains("dcdfile               = ./pro.dcd"));
        assert!(pro.contains("ensemble              = NPT"));
    }

    #[test]
    fn nsteps_must_be_a_multiple_of_each_period() {
        // Dynamics: each period in turn.
        assert_eq!(
            nvt(
                "membrane",
                "min.rst",
                "nvt",
                &Dynamics {
                    rstout_period: 300,
                    ..dynamics()
                }
            )
            .unwrap_err(),
            InpError::NotMultipleOf("rstout_period")
        );
        assert_eq!(
            npt(
                "membrane",
                "min.rst",
                "npt",
                &Dynamics {
                    num_steps: 25000,
                    stoptr_period: 7,
                    ..dynamics()
                }
            )
            .unwrap_err(),
            InpError::NotMultipleOf("stoptr_period")
        );
        assert_eq!(
            nvt(
                "membrane",
                "min.rst",
                "nvt",
                &Dynamics {
                    nbupdate_period: 3,
                    ..dynamics()
                }
            )
            .unwrap_err(),
            InpError::NotMultipleOf("nbupdate_period")
        );

        // Minimization: no stoptr period, but the rest still apply.
        assert_eq!(
            minimization(
                "membrane",
                &Minimization {
                    eneout_period: 7,
                    ..minimization_params()
                }
            )
            .unwrap_err(),
            InpError::NotMultipleOf("eneout_period")
        );
    }

    #[test]
    fn ele_cutoff_scales_for_small_boxes_and_omits_for_large() {
        assert_eq!(ele_cutoff_for_box([14.0, 14.0, 20.0]), None);
        let (rc, rp) = ele_cutoff_for_box([6.887, 6.887, 20.0]).unwrap();
        assert_eq!(rc, 28.0);
        assert_eq!(rp, 33.0);
        assert!(2.0 * rc < 68.87);
        assert!(2.0 * rp < 68.87);
    }

    #[test]
    fn renders_umbrella_restraint_reus_and_mbar_controls() {
        let restr = ZRestraint {
            group1: "rnam:TRP",
            group2: "rnam:POPC",
            function: "DISTMASS",
            direction: "Z2",
            constant: 0.5,
            reference_angstrom: 15.0,
        };
        let min = minimization_with_restraint(
            "../system.top",
            "./system.gro",
            "min",
            &Minimization {
                cg_ele_cutoff: Some((28.0, 33.0)),
                ..minimization_params()
            },
            Some(&restr),
        )
        .unwrap();
        assert!(min.contains("cg_cutoffdist_ele     = 28.0"));
        assert!(min.contains("cg_pairlistdist_ele   = 33.0"));
        assert!(min.contains("group1                = rnam:TRP"));
        assert!(min.contains("group2                = rnam:POPC"));
        assert!(min.contains("function1             = DISTMASS"));
        assert!(min.contains("direction1            = Z2"));
        assert!(min.contains("reference1            = 15.0000"));

        let refs: Vec<f64> = (0..12).map(|i| i as f64 * 2.5).collect();
        let ks = vec![0.5; 12];
        let reus_inp = reus(
            Ensemble::Nvt,
            "./system.top",
            "./win_{}.gro",
            "./min_{}.rst",
            "prod_reus",
            &Dynamics {
                cg_ele_cutoff: Some((28.0, 33.0)),
                ..dynamics()
            },
            &ReusConfig {
                exchange_period: 250,
                iseed: 77,
                group1: "rnam:TRP",
                group2: "rnam:POPC",
                function: "DISTMASS",
                direction: "Z2",
                constants: &ks,
                references_angstrom: &refs,
            },
        )
        .unwrap();
        assert!(reus_inp.contains("type1                 = RESTRAINT"));
        assert!(reus_inp.contains("nreplica1             = 12"));
        assert!(reus_inp.contains("constant1             = \\"));
        assert!(reus_inp.contains("reference1            = \\"));

        let mbar_inp = mbar(&MbarConfig {
            cvfile: "./cv/{}.dat",
            fenefile: "./mbar.fene",
            pmffile: "./mbar.pmf",
            weightfile: Some("./weight_{}.dat"),
            temperature: 323.15,
            nblocks: 5,
            grid_min_angstrom: -1.0,
            grid_max_angstrom: 28.5,
            grid_points: 60,
            constants: &ks,
            references_angstrom: &refs,
        });
        assert!(mbar_inp.contains("nreplica           = 12"));
        assert!(mbar_inp.contains("input_type         = US"));
        assert!(mbar_inp.contains("output_unit        = kcal/mol"));
        assert!(mbar_inp.contains("is_periodic1       = NO"));
        assert!(!mbar_inp.contains("weightfile"));
    }
}
