use crate::inp::boundary::Boundary;
use crate::inp::boundary::BoxSize;
use crate::inp::ensemble::Ensemble;
use crate::inp::output::Output;
use std::fmt;

#[derive(Debug, Default, Clone)]
pub struct InputFile {
    pub input_grotop: String,
    pub input_grocrd: String,
    pub input_rst: Option<String>,
    pub output_rst: Option<Output>,
    pub output_dcd: Option<Output>,
    pub solvent_temperature: f64,
    pub solvent_ionic_strength: f64,
    pub time_step: f64,
    pub num_steps: u64,
    pub output_ene_period: u64,
    pub update_nb_period: u64,
    pub remove_tr_period: u64,
    pub seed: u16,
    pub ensemble: Ensemble,
    pub boundary: Boundary,
}

impl fmt::Display for InputFile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Input
        writeln!(f, "[INPUT]")?;
        writeln!(f, "grotopfile            = {}", self.input_grotop)?;
        writeln!(f, "grocrdfile            = {}", self.input_grocrd)?;
        if let Some(input_rst) = &self.input_rst {
            writeln!(f, "rstfile               = {}", input_rst)?;
        }
        writeln!(f)?;

        // Output
        writeln!(f, "[OUTPUT]")?;
        if let Some(output_rst) = &self.output_rst {
            writeln!(f, "rstfile               = {}", output_rst.path)?;
        }
        if let Some(output_dcd) = &self.output_dcd {
            writeln!(f, "dcdfile               = {}", output_dcd.path)?;
        }
        writeln!(f)?;

        // Energy
        writeln!(f, "[ENERGY]")?;
        writeln!(f, "forcefield            = RESIDCG",)?;
        writeln!(f, "electrostatic         = CUTOFF",)?;
        writeln!(f, "nonb_limiter          = YES",)?;
        writeln!(f, "cg_sol_temperature    = {}", self.solvent_temperature)?;
        writeln!(f, "cg_sol_ionic_strength = {}", self.solvent_ionic_strength)?;
        writeln!(f)?;

        // Dynamics
        writeln!(f, "[DYNAMICS]")?;
        writeln!(f, "integrator            = VVER_CG",)?;
        writeln!(f, "timestep              = {}", self.time_step)?;
        writeln!(f, "nsteps                = {}", self.num_steps)?;
        writeln!(f, "eneout_period         = {}", self.output_ene_period)?;
        if let Some(output_crd) = &self.output_dcd {
            writeln!(f, "crdout_period         = {}", output_crd.period)?;
        }
        if let Some(output_rst) = &self.output_rst {
            writeln!(f, "rstout_period         = {}", output_rst.period)?;
        }
        writeln!(f, "nbupdate_period       = {}", self.update_nb_period)?;
        writeln!(f, "stoptr_period         = {}", self.remove_tr_period)?;
        writeln!(f, "iseed                 = {}", self.seed)?;
        writeln!(f)?;

        // Constraints
        writeln!(f, "[CONSTRAINTS]")?;
        writeln!(f, "rigid_bond            = NO",)?;
        writeln!(f)?;

        // ensemble::Ensemble
        writeln!(f, "[ENSEMBLE]")?;
        match self.ensemble {
            Ensemble::Nve => {
                writeln!(f, "ensemble              = NVE")?;
                writeln!(f, "tpcontrol             = NO")?;
            }
            Ensemble::Nvt {
                temperature,
                gamma_t,
            } => {
                writeln!(f, "ensemble              = NVT")?;
                writeln!(f, "tpcontrol             = LANGEVIN")?;
                writeln!(f, "temperature           = {}", temperature)?;
                writeln!(f, "gamma_t               = {}", gamma_t)?;
            }
            Ensemble::Npt {
                temperature,
                pressure,
                gamma_t,
                gamma_p,
            } => {
                writeln!(f, "ensemble              = NPT")?;
                writeln!(f, "tpcontrol             = LANGEVIN")?;
                writeln!(f, "temperature           = {}", temperature)?;
                writeln!(f, "pressure              = {}", pressure)?;
                writeln!(f, "gamma_t               = {}", gamma_t)?;
                writeln!(f, "gamma_p               = {}", gamma_p)?;
                writeln!(f, "isotropy              = Z-FIXED-SEMI-ISO")?;
            }
        }
        writeln!(f)?;

        // Boundary
        writeln!(f, "[Boundary]",)?;
        match &self.boundary {
            Boundary::NoBc => writeln!(f, "type                  = NOBC")?,
            Boundary::Pbc { box_size } => {
                writeln!(f, "type                  = PBC")?;
                if let Some(BoxSize { x, y, z }) = box_size {
                    writeln!(f, "box_size_x            = {}", x)?;
                    writeln!(f, "box_size_y            = {}", y)?;
                    writeln!(f, "box_size_z            = {}", z)?;
                }
            }
        }

        Ok(())
    }
}
