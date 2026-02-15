use std::fmt;

#[derive(Debug)]
pub struct AtomType {
    pub name: String,
    pub n: u64,
    pub mass: f64,
    pub charge: f64,
    pub ptype: String,
    pub rmin: f64,
    pub eps: f64,
}

impl fmt::Display for AtomType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:>6}{:>4}{:>9.4}{:>9.3}{:>6}{:>9.4}{:>9.4}",
            self.name, self.n, self.mass, self.charge, self.ptype, self.rmin, self.eps
        )?;

        Ok(())
    }
}
