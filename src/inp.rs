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

/// nm → Ångström, the unit GENESIS uses internally and in `box_size`.
const ANGSTROM_PER_NM: f64 = 10.0;

/// Error from generating a control file.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum InpError {
    /// `nsteps` is not an exact multiple of the named output/update period.
    #[error("nsteps must be a multiple of {0}")]
    NotMultipleOf(&'static str),
}

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
}

/// Energy-minimization control file reading the `<structure>.{top,gro}` pair.
///
/// # Errors
///
/// Returns [`InpError::NotMultipleOf`] if `num_steps` is not a multiple of an
/// enabled output or update period.
pub fn minimization(structure: &str, params: &Minimization) -> Result<String, InpError> {
    let Minimization {
        box_size,
        temperature,
        num_steps,
        eneout_period,
        crdout_period,
        rstout_period,
        nbupdate_period,
    } = *params;
    require_multiple(num_steps, eneout_period, "eneout_period")?;
    require_multiple(num_steps, crdout_period, "crdout_period")?;
    require_multiple(num_steps, rstout_period, "rstout_period")?;
    require_multiple(num_steps, nbupdate_period, "nbupdate_period")?;

    let [x, y, z] = box_size.map(|nm| nm * ANGSTROM_PER_NM);
    Ok(format!(
        "\
[INPUT]
grotopfile            = ./{structure}.top
grocrdfile            = ./{structure}.gro

[OUTPUT]
rstfile               = ./min.rst
dcdfile               = ./min.dcd
pdbfile               = ./min.pdb

[ENERGY]
forcefield            = RESIDCG
electrostatic         = CUTOFF
cg_sol_temperature    = {temperature}
cg_sol_ionic_strength = 0.15

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
    ))
}

/// Thermostat ensemble for a dynamics stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ensemble {
    /// Constant volume (vesicles).
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
    dynamics(Ensemble::Nvt, structure, restart_in, basename, params)
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
    dynamics(Ensemble::Npt, structure, restart_in, basename, params)
}

/// Render a dynamics control file over the `<structure>.{top,gro}` pair: it reads
/// `restart_in` and writes `<basename>.{rst,dcd,pdb}`, under the chosen `ensemble`.
fn dynamics(
    ensemble: Ensemble,
    structure: &str,
    restart_in: &str,
    basename: &str,
    params: &Dynamics,
) -> Result<String, InpError> {
    validate_dynamics(params)?;
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
    } = *params;
    let ensemble = ensemble_block(ensemble, temperature);
    Ok(format!(
        "\
[INPUT]
grotopfile            = ./{structure}.top
grocrdfile            = ./{structure}.gro
rstfile               = ./{restart_in}

[OUTPUT]
rstfile               = ./{basename}.rst
dcdfile               = ./{basename}.dcd
pdbfile               = ./{basename}.pdb

[ENERGY]
forcefield            = RESIDCG
electrostatic         = CUTOFF
cg_sol_temperature    = {temperature}
cg_sol_ionic_strength = 0.15

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
    ))
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
}
