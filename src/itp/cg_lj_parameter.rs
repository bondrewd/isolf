use std::fmt;

#[derive(Debug)]
pub struct CgLjParameter {
    pub bead1: String,
    pub bead2: String,
    pub epsilon: f64,
    pub sigma: f64,
    pub cutoff: f64,
}

impl fmt::Display for CgLjParameter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:>6}{:>6}{:>9.3}{:>9.3}{:>9.3}",
            self.bead1, self.bead2, self.epsilon, self.sigma, self.cutoff
        )?;

        Ok(())
    }
}
