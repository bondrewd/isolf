use std::fmt;

#[derive(Debug, Default)]
pub struct ItpFile {
    pub atom_types: Vec<crate::itp::atom_type::AtomType>,
    pub cg_lj_parameters: Vec<crate::itp::cg_lj_parameter::CgLjParameter>,
    pub cg_wca_parameters: Vec<crate::itp::cg_wca_parameter::CgWcaParameter>,
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

        Ok(())
    }
}

impl TryFrom<crate::itp::force_field::ForceField> for ItpFile {
    type Error = crate::error::IsolfError;

    fn try_from(_: crate::itp::force_field::ForceField) -> Result<Self, Self::Error> {
        Ok(crate::itp::itp_file::ItpFile::default())
        //         // Atom types
        //         let mut mass = csv.ff.mass.clone();
        //         mass.sort_by_key(|mass| mass.bead.to_lowercase());

        //         let mut charge = csv.ff.charge.clone();
        //         charge.sort_by_key(|charge| charge.bead.to_lowercase());

        //         if mass.len() != charge.len() {
        //             return Err(IsolfError::invalid_input(
        //                 "atom_type".into(),
        //                 "Mass and charge lists must have the same length".into(),
        //             ));
        //         }

        //         if !mass
        //             .iter()
        //             .zip(charge.iter())
        //             .all(|(mass, charge)| mass.bead == charge.bead)
        //         {
        //             return Err(IsolfError::invalid_input(
        //                 "atom_type".into(),
        //                 "Mass and charge lists must have the same bead type".into(),
        //             ));
        //         }

        //         let atom_types = mass
        //             .iter()
        //             .zip(charge.iter())
        //             .map(|(mass, charge)| IsolfAtomType {
        //                 name: mass.bead.to_uppercase(),
        //                 n: 1,
        //                 mass: mass.m,
        //                 charge: charge.q,
        //                 ptype: "A".into(),
        //                 rmin: 0.0,
        //                 eps: 0.0,
        //             })
        //             .collect();

        //         // CG LJ Parameters
        //         let lj = csv.ff.lj.clone();
        //         let wca = csv.ff.wca.clone();

        //         let mut cg_lj_parameters: Vec<IsolfCgLjParameter> = lj
        //             .iter()
        //             .zip(lj.iter())
        //             .map(|(lj1, lj2)| IsolfCgLjParameter {
        //                 bead1: lj1.bead.to_uppercase(),
        //                 bead2: lj2.bead.to_uppercase(),
        //                 epsilon: (lj1.e * lj2.e).sqrt(),
        //                 sigma: (lj1.s + lj2.s) / 2.0,
        //                 cutoff: 2.5 * (lj1.s + lj2.s) / 2.0,
        //             })
        //             .collect();

        //         let mut lj_wca = lj
        //             .iter()
        //             .flat_map(|lj| {
        //                 wca.iter()
        //                     .filter(|wca| {
        //                         charge
        //                             .iter()
        //                             .any(|charge| charge.bead == wca.bead && charge.q.abs() > 0.0)
        //                     })
        //                     .map(|wca| IsolfCgLjParameter {
        //                         bead1: lj.bead.to_uppercase(),
        //                         bead2: wca.bead.to_uppercase(),
        //                         epsilon: (lj.e * wca.e).sqrt(),
        //                         sigma: (lj.s + wca.s) / 2.0,
        //                         cutoff: 2.5 * (lj.s + wca.s) / 2.0,
        //                     })
        //             })
        //             .collect();

        //         cg_lj_parameters.append(&mut lj_wca);
        //         cg_lj_parameters.sort_by_key(|lj| (lj.bead1.clone(), lj.bead2.clone()));

        //         // CG WCA Parameters
        //         let wca = csv.ff.wca.clone();
        //         let hp = csv.ff.hp.clone();

        //         let mut cg_wca_parameters: Vec<IsolfCgWcaParameter> = wca
        //             .iter()
        //             .flat_map(|wca1| wca.iter().map(move |wca2| (wca1, wca2)))
        //             .filter(|(wca1, wca2)| wca1.bead <= wca2.bead)
        //             .map(|(wca1, wca2)| IsolfCgWcaParameter {
        //                 bead1: wca1.bead.to_uppercase(),
        //                 bead2: wca2.bead.to_uppercase(),
        //                 epsilon: (wca1.e * wca2.e).sqrt(),
        //                 sigma: (wca1.s + wca2.s) / 2.0,
        //             })
        //             .collect();

        //         let mut wca_hp = wca
        //             .iter()
        //             .flat_map(|wca| {
        //                 hp.iter().map(move |hp| IsolfCgWcaParameter {
        //                     bead1: wca.bead.to_uppercase(),
        //                     bead2: hp.bead.to_uppercase(),
        //                     epsilon: (wca.e * hp.e).sqrt(),
        //                     sigma: (wca.s + hp.s) / 2.0,
        //                 })
        //             })
        //             .collect();

        //         cg_wca_parameters.append(&mut wca_hp);
        //         cg_wca_parameters.sort_by_key(|wca| (wca.bead1.clone(), wca.bead2.clone()));

        //         Ok(IsolfItp {
        //             atom_types,
        //             cg_lj_parameters,
        //             cg_wca_parameters,
        //         })
    }
}
