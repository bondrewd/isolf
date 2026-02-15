use crate::error::IsolfError;
use colored::Colorize;
use rand::Rng;

#[derive(Debug, Default, Clone)]
pub struct InputFileBuilder {
    input_grotop: Option<String>,
    input_grocrd: Option<String>,
    input_rst: Option<String>,
    output_rst: Option<crate::inp::output::Output>,
    output_dcd: Option<crate::inp::output::Output>,
    solvent_temperature: Option<f64>,
    solvent_ionic_strength: Option<f64>,
    time_step: Option<f64>,
    num_steps: Option<u64>,
    output_ene_period: Option<u64>,
    update_nb_period: Option<u64>,
    remove_tr_period: Option<u64>,
    seed: Option<u16>,
    ensemble: crate::inp::ensemble::Ensemble,
    boundary: crate::inp::boundary::Boundary,
}

impl InputFileBuilder {
    pub fn build(self) -> Result<crate::inp::inp_file::InputFile, IsolfError> {
        // assert the presence of required fields
        let input_grotop = self
            .input_grotop
            .ok_or_else(|| IsolfError::missing_required_field("input_grotop"))?;

        let input_grocrd = self
            .input_grocrd
            .ok_or_else(|| IsolfError::missing_required_field("input_grocrd"))?;

        let num_steps = self
            .num_steps
            .ok_or_else(|| IsolfError::missing_required_field("num_steps"))?;

        let output_ene_period = self
            .output_ene_period
            .ok_or_else(|| IsolfError::missing_required_field("output_ene_period"))?;

        let update_nb_period = self
            .update_nb_period
            .ok_or_else(|| IsolfError::missing_required_field("update_nb_period"))?;

        let remove_tr_period = self
            .remove_tr_period
            .ok_or_else(|| IsolfError::missing_required_field("remove_tr_period"))?;

        let ensemble = self.ensemble;
        let boundary = self.boundary;

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
            return Err(IsolfError::invalid_field_value(
                "output_ene_period",
                0,
                "period must be greater than zero",
            ));
        }

        // if output_rst is present, the period must be greater than zero
        if let Some(crate::inp::output::Output { period: 0, .. }) = output_rst {
            return Err(IsolfError::invalid_field_value(
                "output_rst",
                0,
                "period must be greater than zero",
            ));
        }

        // if output_dcd is present, the period must be greater than zero
        if let Some(crate::inp::output::Output { period: 0, .. }) = output_dcd {
            return Err(IsolfError::invalid_field_value(
                "output_dcd",
                0,
                "period must be greater than zero",
            ));
        }

        // periods should be divisors of the total number of steps
        if num_steps % output_ene_period != 0 {
            return Err(IsolfError::invalid_field_value(
                "output_ene",
                output_ene_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if let Some(output) = &output_rst
            && num_steps % output.period != 0
        {
            return Err(IsolfError::invalid_field_value(
                "output_rst",
                output.period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if let Some(output) = &output_dcd
            && num_steps % output.period != 0
        {
            return Err(IsolfError::invalid_field_value(
                "output_dcd",
                output.period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if num_steps % update_nb_period != 0 {
            return Err(IsolfError::invalid_field_value(
                "update_nb_period",
                update_nb_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        if num_steps % remove_tr_period != 0 {
            return Err(IsolfError::invalid_field_value(
                "remove_tr_period",
                remove_tr_period,
                "period must be a divisor of the total number of steps",
            ));
        }

        // emit a warning if the solvent temperature is different from
        // the ensemble temperature
        match &ensemble {
            crate::inp::ensemble::Ensemble::Nvt { temperature, .. }
            | crate::inp::ensemble::Ensemble::Npt { temperature, .. } => {
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
        if input_rst.is_some()
            && matches!(
                boundary,
                crate::inp::boundary::Boundary::Pbc { box_size: Some(_) }
            )
        {
            return Err(IsolfError::invalid_field_value(
                "boundary",
                "box_size",
                "box_size is not needed when using a restart file",
            ));
        }

        // when using a PBC boundary, the box size should be present if
        // running a new simulation
        if input_rst.is_none()
            && matches!(
                boundary,
                crate::inp::boundary::Boundary::Pbc { box_size: None }
            )
            && let crate::inp::boundary::Boundary::Pbc { box_size: None } = &boundary
        {
            return Err(IsolfError::invalid_field_value(
                "boundary",
                "box_size",
                "box_size is required when not using a restart file",
            ));
        }

        Ok(crate::inp::inp_file::InputFile {
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

    pub fn output_rst(mut self, path: &str, period: u64) -> Self {
        self.output_rst = Some(crate::inp::output::Output::new(path, period));
        self
    }

    pub fn output_dcd(mut self, path: &str, period: u64) -> Self {
        self.output_dcd = Some(crate::inp::output::Output::new(path, period));
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

    pub fn nvt(mut self, temperature: f64, gamma_t: f64) -> Self {
        self.ensemble = crate::inp::ensemble::Ensemble::nvt(temperature, gamma_t);
        self
    }

    pub fn npt(mut self, temperature: f64, pressure: f64, gamma_t: f64, gamma_p: f64) -> Self {
        self.ensemble =
            crate::inp::ensemble::Ensemble::npt(temperature, pressure, gamma_t, gamma_p);
        self
    }

    pub fn pbc(mut self) -> Self {
        self.boundary = crate::inp::boundary::Boundary::pbc();
        self
    }

    pub fn pbc_with_box_size(mut self, x: f64, y: f64, z: f64) -> Self {
        self.boundary = crate::inp::boundary::Boundary::pbc_with_box_size(x, y, z);
        self
    }
}
