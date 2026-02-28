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
    #[must_use]
    pub const fn nve() -> Self {
        Self::Nve
    }

    #[must_use]
    pub const fn nvt(temperature: f64, gamma_t: f64) -> Self {
        Self::Nvt {
            temperature,
            gamma_t,
        }
    }

    #[must_use]
    pub const fn npt(temperature: f64, pressure: f64, gamma_t: f64, gamma_p: f64) -> Self {
        Self::Npt {
            temperature,
            pressure,
            gamma_t,
            gamma_p,
        }
    }
}
