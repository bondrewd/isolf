use std::fmt;

#[derive(Debug, Default)]
pub struct Bond {
    pub ids: [u64; 2],
    pub k: f64,
    pub r0: f64,
}

impl fmt::Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:>4}{:>4}   1{:>9.4}{:>15.4}",
            self.ids[0], self.ids[1], self.r0, self.k,
        )?;

        Ok(())
    }
}
