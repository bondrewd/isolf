use std::fmt;

#[derive(Debug, Default)]
pub struct Molecule {
    pub molecule_type: crate::itp::molecule_type::MoleculeType,
    pub atoms: Vec<crate::itp::atom::Atom>,
    pub bonds: Vec<crate::itp::bond::Bond>,
    pub angles: Vec<crate::itp::angle::Angle>,
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Molecule type
        writeln!(f, "[ moleculetype ]")?;
        writeln!(f, "; name  nrexcl")?;
        writeln!(f, "{}", self.molecule_type)?;
        writeln!(f)?;

        // Atoms
        writeln!(f, "[ atoms ]")?;
        writeln!(f, "; nr   type  resnr    res   atom   cg   charge     mass")?;
        writeln!(f, ";  -      -      -      -      -    -        e      amu")?;
        for atom in &self.atoms {
            writeln!(f, "{atom}")?;
        }
        writeln!(f)?;

        // Bonds
        writeln!(f, "[ bonds ]")?;
        writeln!(f, ";  i   j   f       eq           coef")?;
        writeln!(f, ";  -   -   -       nm  kJ*nm-2*mol-1")?;
        for bond in &self.bonds {
            writeln!(f, "{bond}")?;
        }
        writeln!(f)?;

        // Angle
        writeln!(f, "[ angles ]")?;
        writeln!(f, ";  i   j   k   f       eq           coef")?;
        writeln!(f, ";  -   -   -   -      deg kJ*rad-2*mol-1")?;
        for angle in &self.angles {
            writeln!(f, "{angle}")?;
        }

        Ok(())
    }
}
