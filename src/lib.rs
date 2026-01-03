use colored::Colorize;
use rand::Rng;
use std::fmt;

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

#[derive(Debug)]
pub enum InputFileBuilderError {
    MissingRequiredField {
        field: String,
    },
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

impl fmt::Display for InputFileBuilderError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InputFileBuilderError::MissingRequiredField { field } => {
                write!(f, "Missing required field '{}'", field.blue().bold())
            }
            InputFileBuilderError::InvalidFieldValue {
                field,
                value,
                reason,
            } => {
                write!(
                    f,
                    "Invalid value '{}' for field '{}': {}",
                    value.blue().bold(),
                    field.blue().bold(),
                    reason.green().bold()
                )
            }
        }
    }
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
        if let Some(output) = &output_rst {
            if output.period == 0 {
                return Err(InputFileBuilderError::invalid_field_value(
                    "output_rst",
                    0,
                    "period must be greater than zero",
                ));
            }
        }

        // if output_dcd is present, the period must be greater than zero
        if let Some(output) = &output_dcd {
            if output.period == 0 {
                return Err(InputFileBuilderError::invalid_field_value(
                    "output_dcd",
                    0,
                    "period must be greater than zero",
                ));
            }
        }

        // periods should be divisors of the total number of steps
        if num_steps % output_ene_period != 0 {
            return Err(InputFileBuilderError::invalid_field_value(
                "output_ene",
                output_ene_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if let Some(output) = &output_rst {
            if num_steps % output.period != 0 {
                return Err(InputFileBuilderError::invalid_field_value(
                    "output_rst",
                    output.period,
                    "period must be a divisor of the total number of steps",
                ));
            }
        }

        if let Some(output) = &output_dcd {
            if num_steps % output.period != 0 {
                return Err(InputFileBuilderError::invalid_field_value(
                    "output_dcd",
                    output.period,
                    "period must be a divisor of the total number of steps",
                ));
            }
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
        if input_rst.is_none() && matches!(boundary, Boundary::Pbc { box_size: None }) {
            if let Boundary::Pbc { box_size: None } = &boundary {
                return Err(InputFileBuilderError::invalid_field_value(
                    "boundary",
                    "box_size",
                    "box_size is required when not using a restart file",
                ));
            }
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
        write!(f, "[INPUT]\n",)?;
        write!(f, "grotopfile            = {}\n", self.input_grotop)?;
        write!(f, "grocrdfile            = {}\n", self.input_grocrd)?;
        if let Some(input_rst) = &self.input_rst {
            write!(f, "rstfile               = {}\n", input_rst)?;
        }
        write!(f, "\n",)?;

        // Output
        write!(f, "[OUTPUT]\n",)?;
        if let Some(output_rst) = &self.output_rst {
            write!(f, "rstfile               = {}\n", output_rst.path)?;
        }
        if let Some(output_dcd) = &self.output_dcd {
            write!(f, "dcdfile               = {}\n", output_dcd.path)?;
        }
        write!(f, "\n",)?;

        // Energy
        write!(f, "[ENERGY]\n",)?;
        write!(f, "forcefield            = RESIDCG\n",)?;
        write!(f, "electrostatic         = CUTOFF\n",)?;
        write!(f, "nonb_limiter          = YES\n",)?;
        write!(f, "cg_sol_temperature    = {}\n", self.solvent_temperature)?;
        write!(
            f,
            "cg_sol_ionic_strength = {}\n",
            self.solvent_ionic_strength
        )?;
        write!(f, "\n",)?;

        // Dynamics
        write!(f, "[DYNAMICS]\n",)?;
        write!(f, "integrator            = VVER_CG\n",)?;
        write!(f, "timestep              = {}\n", self.time_step)?;
        write!(f, "nsteps                = {}\n", self.num_steps)?;
        write!(f, "eneout_period         = {}\n", self.output_ene_period)?;
        if let Some(output_crd) = &self.output_dcd {
            write!(f, "crdout_period         = {}\n", output_crd.period)?;
        }
        if let Some(output_rst) = &self.output_rst {
            write!(f, "rstout_period         = {}\n", output_rst.period)?;
        }
        write!(f, "nbupdate_period       = {}\n", self.update_nb_period)?;
        write!(f, "stoptr_period         = {}\n", self.remove_tr_period)?;
        write!(f, "iseed                 = {}\n", self.seed)?;
        write!(f, "\n",)?;

        // Constraints
        write!(f, "[CONSTRAINTS]\n",)?;
        write!(f, "rigid_bond            = NO\n",)?;
        write!(f, "\n",)?;

        // Ensemble
        write!(f, "[ENSEMBLE]\n",)?;
        match self.ensemble {
            Ensemble::Nve => {
                write!(f, "ensemble              = NVE\n")?;
                write!(f, "tpcontrol             = NO\n")?;
            }
            Ensemble::Nvt {
                temperature,
                gamma_t,
            } => {
                write!(f, "ensemble              = NVT\n")?;
                write!(f, "tpcontrol             = LANGEVIN\n")?;
                write!(f, "temperature           = {}\n", temperature)?;
                write!(f, "gamma_t               = {}\n", gamma_t)?;
            }
            Ensemble::Npt {
                temperature,
                pressure,
                gamma_t,
                gamma_p,
            } => {
                write!(f, "ensemble              = NPT\n")?;
                write!(f, "tpcontrol             = LANGEVIN\n")?;
                write!(f, "temperature           = {}\n", temperature)?;
                write!(f, "pressure              = {}\n", pressure)?;
                write!(f, "gamma_t               = {}\n", gamma_t)?;
                write!(f, "gamma_p               = {}\n", gamma_p)?;
                write!(f, "isotropy              = Z-FIXED-SEMI-ISO\n")?;
            }
        }
        write!(f, "\n",)?;

        // Boundary
        write!(f, "[Boundary]\n",)?;
        match &self.boundary {
            Boundary::NoBc => write!(f, "type                  = NOBC\n")?,
            Boundary::Pbc { box_size } => {
                write!(f, "type                  = PBC\n")?;
                if let Some(box_size) = box_size {
                    write!(f, "box_size_x            = {}\n", box_size.x)?;
                    write!(f, "box_size_y            = {}\n", box_size.y)?;
                    write!(f, "box_size_z            = {}\n", box_size.z)?;
                }
            }
        }

        Ok(())
    }
}
