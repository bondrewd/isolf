use std::fmt;

#[derive(Debug, Default)]
pub struct Angle {
    pub ids: [u64; 3],
    pub k: f64,
    pub t0: f64,
}

impl fmt::Display for Angle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:>4}{:>4}{:>4}   1{:>9.4}{:>15.4}",
            self.ids[0], self.ids[1], self.ids[2], self.t0, self.k,
        )?;

        Ok(())
    }
}
