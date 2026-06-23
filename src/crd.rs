//! CHARMM coordinate (`.crd`) rendering of a particle list.
//!
//! Coordinates are written in ångström (nm × 10). The standard fixed-column
//! format is used until the system exceeds its five-digit / four-character
//! fields, at which point the extended format is emitted automatically. A
//! [`Membrane`] renders through the [`render`] core.

use std::fmt::Write;

use crate::membrane::{Membrane, Particle};

/// nm → ångström conversion for CHARMM coordinates.
const ANGSTROM_PER_NM: f64 = 10.0;

/// Segment name used for every lipid molecule.
const SEGMENT: &str = "MEMB";

/// Largest atom count the standard (five-digit) format can represent.
const MAX_STANDARD_ATOMS: usize = 99_999;

/// Largest residue id the standard (four-character) resid field can represent.
const MAX_STANDARD_RESIDUES: usize = 9_999;

/// Render a particle list as a CHARMM `.crd` coordinate file. `num_residues` is
/// the molecule count, used only to pick the standard or extended layout.
pub fn render(particles: &[Particle], name: &str, num_residues: usize) -> String {
    let extended = particles.len() > MAX_STANDARD_ATOMS || num_residues > MAX_STANDARD_RESIDUES;

    let mut s = String::new();
    let _ = writeln!(s, "* {name}");
    let _ = writeln!(s, "*");
    if extended {
        let _ = writeln!(s, "{:>10}  EXT", particles.len());
    } else {
        let _ = writeln!(s, "{:>5}", particles.len());
    }

    for particle in particles {
        let [x, y, z] = particle.position.map(|c| c * ANGSTROM_PER_NM);
        if extended {
            let _ = writeln!(
                s,
                "{:>10}{:>10}  {:<8}  {:<8}{:>20.10}{:>20.10}{:>20.10}  {:<8}  {:<8}{:>20.10}",
                particle.atom_id,
                particle.molecule_id,
                particle.residue,
                particle.bead,
                x,
                y,
                z,
                SEGMENT,
                particle.molecule_id,
                0.0,
            );
        } else {
            let _ = writeln!(
                s,
                "{:>5}{:>5} {:<4} {:<4}{:>10.5}{:>10.5}{:>10.5} {:<4} {:<4}{:>10.5}",
                particle.atom_id,
                particle.molecule_id,
                particle.residue,
                particle.bead,
                x,
                y,
                z,
                SEGMENT,
                particle.molecule_id,
                0.0,
            );
        }
    }
    s
}

impl Membrane {
    /// Render this membrane as a CHARMM `.crd` coordinate file.
    pub fn to_crd(&self) -> String {
        render(&self.particles, &self.name, self.total_lipids())
    }
}

#[cfg(test)]
mod tests {
    use crate::composition::Composition;
    use crate::force_field::ForceField;
    use crate::membrane::{BuildOptions, DEFAULT_NAME, Membrane, Sizing};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn membrane(count: usize) -> Membrane {
        let force_field = ForceField::isolf();
        let leaflet = Composition::parse("DLPC=1").unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        Membrane::build(
            &force_field,
            DEFAULT_NAME,
            &leaflet,
            &leaflet,
            Sizing::Count {
                upper: count,
                lower: count,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap()
    }

    #[test]
    fn standard_header_and_atom_count() {
        let membrane = membrane(4);
        let crd = membrane.to_crd();
        let lines: Vec<&str> = crd.lines().collect();

        assert_eq!(lines[0], "* CG membrane model");
        assert_eq!(lines[1], "*");
        // 8 DLPC × 5 beads = 40 atoms, standard 5-wide count.
        assert_eq!(lines[2].trim(), "40");
        assert_eq!(lines.iter().filter(|l| l.contains("DLPC")).count(), 40);
    }

    #[test]
    fn coordinates_are_in_angstrom() {
        let membrane = membrane(4);
        let crd = membrane.to_crd();
        // A DLPC head bead reaches a few nm in z, i.e. tens of ångström. Confirm
        // a coordinate magnitude consistent with ångström, not nm.
        let max_coord = membrane
            .particles
            .iter()
            .flat_map(|p| p.position)
            .fold(0.0_f64, |m, c| m.max(c.abs()));
        assert!(crd.contains(&format!("{:.5}", max_coord * 10.0)));
    }
}
