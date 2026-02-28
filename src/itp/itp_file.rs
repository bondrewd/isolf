use crate::error::IsolfError;
use crate::itp::angle::Angle;
use crate::itp::atom::Atom;
use crate::itp::atom_type::AtomType;
use crate::itp::bond::Bond;
use crate::itp::cg_lj_parameter::CgLjParameter;
use crate::itp::cg_wca_parameter::CgWcaParameter;
use crate::itp::force_field::Bead;
use crate::itp::force_field::ForceField;
use crate::itp::molecule::Molecule;
use std::collections::HashMap;
use std::fmt;

type AtomTypes = Vec<AtomType>;
type CgLjParameters = Vec<CgLjParameter>;
type CgWcaParameters = Vec<CgWcaParameter>;
type Molecules = Vec<Molecule>;

impl From<ForceField> for AtomTypes {
    fn from(ff: ForceField) -> Self {
        ff.atoms
            .iter()
            .map(|atom| AtomType {
                name: atom.name.to_uppercase(),
                n: 1,
                mass: atom.m,
                charge: atom.q,
                ptype: "A".into(),
                rmin: 0.0,
                eps: 0.0,
            })
            .collect()
    }
}

impl From<ForceField> for CgLjParameters {
    fn from(ff: ForceField) -> Self {
        let polar_atoms: Vec<&Bead> = ff.atoms.iter().filter(|atom| atom.p).collect();
        let charged_atoms: Vec<&Bead> = ff.atoms.iter().filter(|atom| atom.q != 0.0).collect();
        let cg_lj_parameters_polar_polar = polar_atoms.iter().flat_map(|atom1| {
            polar_atoms.iter().map(|atom2| CgLjParameter {
                bead1: atom1.name.to_uppercase(),
                bead2: atom2.name.to_uppercase(),
                epsilon: (atom1.e * atom2.e).sqrt(),
                sigma: f64::midpoint(atom1.s, atom2.s),
                cutoff: 2.5 * (atom1.s + atom2.s) / 2.0,
            })
        });
        let cg_lj_parameters_polar_charged: Self = polar_atoms
            .iter()
            .flat_map(|atom1| {
                charged_atoms.iter().map(|atom2| CgLjParameter {
                    bead1: atom1.name.to_uppercase(),
                    bead2: atom2.name.to_uppercase(),
                    epsilon: (atom1.e * atom2.e).sqrt(),
                    sigma: f64::midpoint(atom1.s, atom2.s),
                    cutoff: 2.5 * (atom1.s + atom2.s) / 2.0,
                })
            })
            .collect();
        cg_lj_parameters_polar_polar
            .chain(cg_lj_parameters_polar_charged)
            .collect()
    }
}

impl From<ForceField> for CgWcaParameters {
    fn from(ff: ForceField) -> Self {
        let polar_atoms: Vec<&Bead> = ff.atoms.iter().filter(|atom| atom.p).collect();
        let no_tail_atoms: Vec<&Bead> = ff.atoms.iter().filter(|atom| atom.w.is_none()).collect();
        let tail_atoms: Vec<&Bead> = ff.atoms.iter().filter(|atom| atom.w.is_some()).collect();
        let cg_wca_parameters_no_tail_no_tail = no_tail_atoms.iter().flat_map(|atom1| {
            no_tail_atoms.iter().map(|atom2| CgWcaParameter {
                bead1: atom1.name.to_uppercase(),
                bead2: atom2.name.to_uppercase(),
                epsilon: (atom1.e * atom2.e).sqrt(),
                sigma: f64::midpoint(atom1.s, atom2.s),
            })
        });
        let cg_wca_parameters_no_tail_tail: Self = polar_atoms
            .iter()
            .flat_map(|atom1| {
                tail_atoms.iter().map(|atom2| CgWcaParameter {
                    bead1: atom1.name.to_uppercase(),
                    bead2: atom2.name.to_uppercase(),
                    epsilon: (atom1.e * atom2.e).sqrt(),
                    sigma: f64::midpoint(atom1.s, atom2.s),
                })
            })
            .collect();
        cg_wca_parameters_no_tail_no_tail
            .chain(cg_wca_parameters_no_tail_tail)
            .collect()
    }
}

impl From<ForceField> for Molecules {
    fn from(ff: ForceField) -> Self {
        let atoms_map: HashMap<String, (f64, f64)> = ff
            .clone()
            .atoms
            .into_iter()
            .map(|atom| (atom.name, (atom.m, atom.q)))
            .collect();
        ff.lipids
            .iter()
            .map(|lipid| {
                let mut molecule = Molecule::default();
                // Molecule type
                molecule.molecule_type.name = lipid.name.to_uppercase();
                // Molecule atoms
                molecule.atoms = lipid
                    .atoms
                    .iter()
                    .enumerate()
                    .map(|(id, atom)| {
                        let (mass, charge) = atoms_map.get(&atom.name).unwrap();
                        Atom {
                            atom_name: atom.name.to_uppercase(),
                            molecule_name: lipid.name.to_uppercase(),
                            q: *charge,
                            m: *mass,
                            id: 1 + id as u64,
                        }
                    })
                    .collect();
                // Molecule bonds
                molecule.bonds = lipid
                    .bonds
                    .iter()
                    .map(|bond| Bond {
                        ids: [bond.ids[0] + 1, bond.ids[1] + 1],
                        k: bond.k,
                        r0: bond.r0,
                    })
                    .collect();
                // Molecule angles
                molecule.angles = lipid
                    .angles
                    .iter()
                    .map(|angle| Angle {
                        ids: [angle.ids[0] + 1, angle.ids[1] + 1, angle.ids[2] + 1],
                        k: angle.k,
                        t0: angle.t0,
                    })
                    .collect();

                molecule
            })
            .collect()
    }
}

#[derive(Debug, Default)]
pub struct ItpFile {
    pub atom_types: AtomTypes,
    pub cg_lj_parameters: CgLjParameters,
    pub cg_wca_parameters: CgWcaParameters,
    pub molecules: Molecules,
}

impl fmt::Display for ItpFile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Atom types
        writeln!(f, "[ atomtypes ]")?;
        writeln!(f, "; name   n     mass   charge ptype     rmin      eps")?;
        writeln!(f, ";    -   -    g/mol        e     -       nm   kJ/mol")?;
        for atom_type in &self.atom_types {
            writeln!(f, "{atom_type}")?;
        }
        writeln!(f)?;

        // CG LJ Parameters
        writeln!(f, "[ cg_LJ_parameters ]")?;
        writeln!(f, "; name  name  epsilon    sigma  cut-off")?;
        writeln!(f, ";    -     -   kJ/mol       nm       nm")?;
        for cg_lj_parameter in &self.cg_lj_parameters {
            writeln!(f, "{cg_lj_parameter}")?;
        }
        writeln!(f)?;

        // CG WCA Parameters
        writeln!(f, "[ cg_WCA_parameters ]")?;
        writeln!(f, "; name  name  epsilon    sigma")?;
        writeln!(f, ";    -     -   kJ/mol       nm")?;
        for cg_wca_parameter in &self.cg_wca_parameters {
            writeln!(f, "{cg_wca_parameter}")?;
        }

        // Molecules
        for molecule in &self.molecules {
            writeln!(f, "{molecule}")?;
        }

        Ok(())
    }
}

impl TryFrom<ForceField> for ItpFile {
    type Error = IsolfError;

    fn try_from(ff: ForceField) -> Result<Self, Self::Error> {
        let itp_file = Self {
            atom_types: AtomTypes::from(ff.clone()),
            cg_lj_parameters: CgLjParameters::from(ff.clone()),
            cg_wca_parameters: CgWcaParameters::from(ff.clone()),
            molecules: Molecules::from(ff),
        };

        Ok(itp_file)
    }
}
