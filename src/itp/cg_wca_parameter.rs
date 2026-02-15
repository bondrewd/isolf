use std::fmt;

#[derive(Debug)]
pub struct CgWcaParameter {
    pub bead1: String,
    pub bead2: String,
    pub epsilon: f64,
    pub sigma: f64,
}

impl fmt::Display for CgWcaParameter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:>6}{:>6}{:>9.3}{:>9.3}",
            self.bead1, self.bead2, self.epsilon, self.sigma
        )?;

        Ok(())
    }
}
