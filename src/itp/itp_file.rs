use std::fmt;

#[derive(Debug, Default)]
pub struct ItpFile {
    pub atom_types: Vec<crate::itp::atom_type::AtomType>,
    pub cg_lj_parameters: Vec<crate::itp::cg_lj_parameter::CgLjParameter>,
    pub cg_wca_parameters: Vec<crate::itp::cg_wca_parameter::CgWcaParameter>,
    pub molecules: Vec<crate::itp::molecule::Molecule>,
}

impl fmt::Display for ItpFile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Atom types
        writeln!(f, "[ atomtypes ]")?;
        writeln!(f, "; name   n     mass   charge ptype     rmin      eps")?;
        writeln!(f, ";    -   -    g/mol        e     -       nm   kJ/mol")?;
        for atom_type in &self.atom_types {
            writeln!(f, "{}", atom_type)?;
        }
        writeln!(f)?;

        // CG LJ Parameters
        writeln!(f, "[ cg_LJ_parameters ]")?;
        writeln!(f, "; name  name  epsilon    sigma  cut-off")?;
        writeln!(f, ";    -     -   kJ/mol       nm       nm")?;
        for cg_lj_parameter in &self.cg_lj_parameters {
            writeln!(f, "{}", cg_lj_parameter)?;
        }
        writeln!(f)?;

        // CG WCA Parameters
        writeln!(f, "[ cg_WCA_parameters ]")?;
        writeln!(f, "; name  name  epsilon    sigma")?;
        writeln!(f, ";    -     -   kJ/mol       nm")?;
        for cg_wca_parameter in &self.cg_wca_parameters {
            writeln!(f, "{}", cg_wca_parameter)?;
        }

        // Molecules
        for molecule in &self.molecules {
            writeln!(f, "{}", molecule)?;
            writeln!(f)?;
        }

        Ok(())
    }
}

impl TryFrom<crate::itp::force_field::ForceField> for ItpFile {
    type Error = crate::error::IsolfError;

    fn try_from(ff: crate::itp::force_field::ForceField) -> Result<Self, Self::Error> {
        // Atom types
        let atom_types = ff
            .atoms
            .iter()
            .map(|atom| crate::itp::atom_type::AtomType {
                name: atom.name.to_uppercase(),
                n: 1,
                mass: atom.m,
                charge: atom.q,
                ptype: "A".into(),
                rmin: 0.0,
                eps: 0.0,
            })
            .collect();

        // CG LJ Parameters
        let polar_atoms: Vec<&crate::itp::force_field::Atom> =
            ff.atoms.iter().filter(|atom| atom.p).collect();
        let charged_atoms: Vec<&crate::itp::force_field::Atom> =
            ff.atoms.iter().filter(|atom| atom.q != 0.0).collect();
        let cg_lj_parameters_polar_polar: Vec<crate::itp::cg_lj_parameter::CgLjParameter> =
            polar_atoms
                .iter()
                .flat_map(|atom1| {
                    polar_atoms
                        .iter()
                        .map(|atom2| crate::itp::cg_lj_parameter::CgLjParameter {
                            bead1: atom1.name.to_uppercase(),
                            bead2: atom2.name.to_uppercase(),
                            epsilon: (atom1.e * atom2.e).sqrt(),
                            sigma: (atom1.s + atom2.s) / 2.0,
                            cutoff: 2.5 * (atom1.s + atom2.s) / 2.0,
                        })
                })
                .collect();
        let cg_lj_parameters_polar_charged: Vec<crate::itp::cg_lj_parameter::CgLjParameter> =
            polar_atoms
                .iter()
                .flat_map(|atom1| {
                    charged_atoms
                        .iter()
                        .map(|atom2| crate::itp::cg_lj_parameter::CgLjParameter {
                            bead1: atom1.name.to_uppercase(),
                            bead2: atom2.name.to_uppercase(),
                            epsilon: (atom1.e * atom2.e).sqrt(),
                            sigma: (atom1.s + atom2.s) / 2.0,
                            cutoff: 2.5 * (atom1.s + atom2.s) / 2.0,
                        })
                })
                .collect();
        let cg_lj_parameters = cg_lj_parameters_polar_polar
            .into_iter()
            .chain(cg_lj_parameters_polar_charged.into_iter())
            .collect();

        // CG WCA Parameters
        let no_tail_atoms: Vec<&crate::itp::force_field::Atom> =
            ff.atoms.iter().filter(|atom| atom.w.is_none()).collect();
        let tail_atoms: Vec<&crate::itp::force_field::Atom> =
            ff.atoms.iter().filter(|atom| atom.w.is_some()).collect();
        let cg_wca_parameters_no_tail_no_tail: Vec<crate::itp::cg_wca_parameter::CgWcaParameter> =
            no_tail_atoms
                .iter()
                .flat_map(|atom1| {
                    no_tail_atoms
                        .iter()
                        .map(|atom2| crate::itp::cg_wca_parameter::CgWcaParameter {
                            bead1: atom1.name.to_uppercase(),
                            bead2: atom2.name.to_uppercase(),
                            epsilon: (atom1.e * atom2.e).sqrt(),
                            sigma: (atom1.s + atom2.s) / 2.0,
                        })
                })
                .collect();
        let cg_wca_parameters_no_tail_tail: Vec<crate::itp::cg_wca_parameter::CgWcaParameter> =
            polar_atoms
                .iter()
                .flat_map(|atom1| {
                    tail_atoms
                        .iter()
                        .map(|atom2| crate::itp::cg_wca_parameter::CgWcaParameter {
                            bead1: atom1.name.to_uppercase(),
                            bead2: atom2.name.to_uppercase(),
                            epsilon: (atom1.e * atom2.e).sqrt(),
                            sigma: (atom1.s + atom2.s) / 2.0,
                        })
                })
                .collect();
        let cg_wca_parameters = cg_wca_parameters_no_tail_no_tail
            .into_iter()
            .chain(cg_wca_parameters_no_tail_tail.into_iter())
            .collect();

        let mut itp_file = crate::itp::itp_file::ItpFile::default();
        itp_file.atom_types = atom_types;
        itp_file.cg_lj_parameters = cg_lj_parameters;
        itp_file.cg_wca_parameters = cg_wca_parameters;

        Ok(itp_file)
    }
}
