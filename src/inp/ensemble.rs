#[derive(Debug, Default, Clone)]
pub enum Ensemble {
    #[default]
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
