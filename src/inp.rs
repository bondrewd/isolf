use colored::Colorize;
use rand::Rng;
use std::fmt;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum InputFileBuilderError {
    #[error("Missing required field '{field}'")]
    MissingRequiredField { field: String },

    #[error("Invalid value '{value}' for field '{field}': {reason}")]
    InvalidFieldValue {
        field: String,
        value: String,
        reason: String,
    },
}

impl InputFileBuilderError {
    pub fn missing_required_field<F: ToString>(field: F) -> Self {
        InputFileBuilderError::MissingRequiredField {
            field: field.to_string(),
        }
    }

    pub fn invalid_field_value<F: ToString, V: ToString, R: ToString>(
        field: F,
        value: V,
        reason: R,
    ) -> Self {
        InputFileBuilderError::InvalidFieldValue {
            field: field.to_string(),
            value: value.to_string(),
            reason: reason.to_string(),
        }
    }
}

#[derive(Debug)]
pub struct Output {
    path: String,
    period: u64,
}

impl Output {
    pub fn new(path: &str, period: u64) -> Self {
        Output {
            path: path.to_string(),
            period,
        }
    }
}

#[derive(Debug)]
pub enum Ensemble {
    Nve,
    Nvt {
        temperature: f64,
        gamma_t: f64,
    },
    Npt {
        temperature: f64,
        pressure: f64,
        gamma_t: f64,
        gamma_p: f64,
    },
}

impl Ensemble {
    pub fn nve() -> Self {
        Ensemble::Nve
    }

    pub fn nvt(temperature: f64, gamma_t: f64) -> Self {
        Ensemble::Nvt {
            temperature,
            gamma_t,
        }
    }

    pub fn npt(temperature: f64, pressure: f64, gamma_t: f64, gamma_p: f64) -> Self {
        Ensemble::Npt {
            temperature,
            pressure,
            gamma_t,
            gamma_p,
        }
    }
}

#[derive(Debug)]
pub struct BoxSize {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl BoxSize {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        BoxSize { x, y, z }
    }
}

#[derive(Debug)]
pub enum Boundary {
    NoBc,
    Pbc { box_size: Option<BoxSize> },
}

impl Boundary {
    pub fn nobc() -> Self {
        Boundary::NoBc
    }

    pub fn pbc() -> Self {
        Boundary::Pbc { box_size: None }
    }

    pub fn pbc_with_box_size(x: f64, y: f64, z: f64) -> Self {
        Boundary::Pbc {
            box_size: Some(BoxSize::new(x, y, z)),
        }
    }
}

#[derive(Debug, Default)]
pub struct InputFileBuilder {
    input_grotop: Option<String>,
    input_grocrd: Option<String>,
    input_rst: Option<String>,
    output_rst: Option<Output>,
    output_dcd: Option<Output>,
    solvent_temperature: Option<f64>,
    solvent_ionic_strength: Option<f64>,
    time_step: Option<f64>,
    num_steps: Option<u64>,
    output_ene_period: Option<u64>,
    update_nb_period: Option<u64>,
    remove_tr_period: Option<u64>,
    seed: Option<u16>,
    ensemble: Option<Ensemble>,
    boundary: Option<Boundary>,
}

impl InputFileBuilder {
    pub fn build(self) -> Result<InputFile, InputFileBuilderError> {
        // assert the presence of required fields
        let input_grotop = self
            .input_grotop
            .ok_or_else(|| InputFileBuilderError::missing_required_field("input_grotop"))?;

        let input_grocrd = self
            .input_grocrd
            .ok_or_else(|| InputFileBuilderError::missing_required_field("input_grocrd"))?;

        let num_steps = self
            .num_steps
            .ok_or_else(|| InputFileBuilderError::missing_required_field("num_steps"))?;

        let output_ene_period = self
            .output_ene_period
            .ok_or_else(|| InputFileBuilderError::missing_required_field("output_ene_period"))?;

        let update_nb_period = self
            .update_nb_period
            .ok_or_else(|| InputFileBuilderError::missing_required_field("update_nb_period"))?;

        let remove_tr_period = self
            .remove_tr_period
            .ok_or_else(|| InputFileBuilderError::missing_required_field("remove_tr_period"))?;

        let ensemble = self
            .ensemble
            .ok_or_else(|| InputFileBuilderError::missing_required_field("ensemble"))?;

        let boundary = self
            .boundary
            .ok_or_else(|| InputFileBuilderError::missing_required_field("boundary"))?;

        // set default values for optional fields
        let time_step = self.time_step.unwrap_or(0.01);
        let solvent_temperature = self.solvent_temperature.unwrap_or(303.15);
        let solvent_ionic_strength = self.solvent_ionic_strength.unwrap_or(0.15);
        let seed = self.seed.unwrap_or_else(|| rand::rng().random());

        // get optional fields
        let input_rst = self.input_rst;
        let output_rst = self.output_rst;
        let output_dcd = self.output_dcd;

        // output_ene_period needs to be greater than zero
        if output_ene_period == 0 {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_ene_period",
                0,
                "period must be greater than zero",
            ));
        }

        // if output_rst is present, the period must be greater than zero
        if let Some(Output { period: 0, .. }) = output_rst {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_rst",
                0,
                "period must be greater than zero",
            ));
        }

        // if output_dcd is present, the period must be greater than zero
        if let Some(Output { period: 0, .. }) = output_dcd {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_dcd",
                0,
                "period must be greater than zero",
            ));
        }

        // periods should be divisors of the total number of steps
        if num_steps % output_ene_period != 0 {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_ene",
                output_ene_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if let Some(output) = &output_rst
            && num_steps % output.period != 0
        {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_rst",
                output.period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if let Some(output) = &output_dcd
            && num_steps % output.period != 0
        {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_dcd",
                output.period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if num_steps % update_nb_period != 0 {
            return Err(InputFileBuilderError::invalid_field_value(
                "update_nb_period",
                update_nb_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if num_steps % remove_tr_period != 0 {
            return Err(InputFileBuilderError::invalid_field_value(
                "remove_tr_period",
                remove_tr_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        // emit a warning if the solvent temperature is different from
        // the ensemble temperature
        match &ensemble {
            Ensemble::Nvt { temperature, .. } | Ensemble::Npt { temperature, .. } => {
                if *temperature != solvent_temperature {
                    eprintln!(
                        "{} solvent temperature ({}) is different from ensemble temperature ({})",
                        "Warning:".yellow().bold(),
                        solvent_temperature,
                        temperature
                    );
                }
            }
            _ => (),
        }

        // when using a PBC boundary, the box size should not be present
        // if continuing a simulation from a restart file
        if input_rst.is_some() && matches!(boundary, Boundary::Pbc { box_size: Some(_) }) {
            return Err(InputFileBuilderError::invalid_field_value(
                "boundary",
                "box_size",
                "box_size is not needed when using a restart file",
            ));
        }

        // when using a PBC boundary, the box size should be present if
        // running a new simulation
        if input_rst.is_none()
            && matches!(boundary, Boundary::Pbc { box_size: None })
            && let Boundary::Pbc { box_size: None } = &boundary
        {
            return Err(InputFileBuilderError::invalid_field_value(
                "boundary",
                "box_size",
                "box_size is required when not using a restart file",
            ));
        }

        Ok(InputFile {
            input_grotop,
            input_grocrd,
            input_rst,
            output_rst,
            output_dcd,
            solvent_temperature,
            solvent_ionic_strength,
            time_step,
            num_steps,
            output_ene_period,
            update_nb_period,
            remove_tr_period,
            seed,
            ensemble,
            boundary,
        })
    }

    pub fn input_grotop(mut self, input_grotop: &str) -> Self {
        self.input_grotop = Some(input_grotop.into());
        self
    }

    pub fn input_grocrd(mut self, input_grocrd: &str) -> Self {
        self.input_grocrd = Some(input_grocrd.into());
        self
    }

    pub fn input_rst(mut self, input_rst: &str) -> Self {
        self.input_rst = Some(input_rst.into());
        self
    }

    pub fn output_rst(mut self, output_rst: Output) -> Self {
        self.output_rst = Some(output_rst);
        self
    }

    pub fn output_dcd(mut self, output_dcd: Output) -> Self {
        self.output_dcd = Some(output_dcd);
        self
    }

    pub fn solvent_temperature(mut self, solvent_temperature: f64) -> Self {
        self.solvent_temperature = Some(solvent_temperature);
        self
    }

    pub fn solvent_ionic_strength(mut self, solvent_ionic_strength: f64) -> Self {
        self.solvent_ionic_strength = Some(solvent_ionic_strength);
        self
    }

    pub fn time_step(mut self, time_step: f64) -> Self {
        self.time_step = Some(time_step);
        self
    }

    pub fn num_steps(mut self, num_steps: u64) -> Self {
        self.num_steps = Some(num_steps);
        self
    }

    pub fn output_ene_period(mut self, output_ene_period: u64) -> Self {
        self.output_ene_period = Some(output_ene_period);
        self
    }

    pub fn update_nb_period(mut self, update_nb_period: u64) -> Self {
        self.update_nb_period = Some(update_nb_period);
        self
    }

    pub fn remove_tr_period(mut self, remove_tr_period: u64) -> Self {
        self.remove_tr_period = Some(remove_tr_period);
        self
    }

    pub fn seed(mut self, seed: u16) -> Self {
        self.seed = Some(seed);
        self
    }

    pub fn ensemble(mut self, ensemble: Ensemble) -> Self {
        self.ensemble = Some(ensemble);
        self
    }

    pub fn boundary(mut self, boundary: Boundary) -> Self {
        self.boundary = Some(boundary);
        self
    }
}

pub struct InputFile {
    input_grotop: String,
    input_grocrd: String,
    input_rst: Option<String>,
    output_rst: Option<Output>,
    output_dcd: Option<Output>,
    solvent_temperature: f64,
    solvent_ionic_strength: f64,
    time_step: f64,
    num_steps: u64,
    output_ene_period: u64,
    update_nb_period: u64,
    remove_tr_period: u64,
    seed: u16,
    ensemble: Ensemble,
    boundary: Boundary,
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

        // Ensemble
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
                if let Some(box_size) = box_size {
                    writeln!(f, "box_size_x            = {}", box_size.x)?;
                    writeln!(f, "box_size_y            = {}", box_size.y)?;
                    writeln!(f, "box_size_z            = {}", box_size.z)?;
                }
            }
        }

        Ok(())
    }
}
