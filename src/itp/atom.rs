use std::fmt;

#[derive(Debug, Default)]
pub struct Atom {
    pub atom_name: String,
    pub molecule_name: String,
    pub q: f64,
    pub m: f64,
    pub id: u64,
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:>4}{:>7}      1{:>7}{:>7}    1{:>9.3}{:>9.4}",
            self.id, self.atom_name, self.molecule_name, self.atom_name, self.q, self.m,
        )?;

        Ok(())
    }
}
