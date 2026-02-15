use std::fmt;

#[derive(Debug, Default)]
pub struct MoleculeType {
    pub name: String,
}

impl fmt::Display for MoleculeType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:>6}              2", self.name,)?;

        Ok(())
    }
}
