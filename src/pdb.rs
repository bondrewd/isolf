//! PDB coordinate (`.pdb`) rendering of a particle list.
//!
//! Particles are written in order with their molecule and atom numbers.
//! Coordinates are converted from nm to ångström (×10), as the PDB format
//! expects. A [`Membrane`] renders through the [`render`] core.

use std::fmt::Write;

use crate::membrane::{Membrane, Particle};

/// nm → ångström conversion for PDB coordinates.
const ANGSTROM_PER_NM: f64 = 10.0;

/// Render a particle list as a `.pdb` coordinate file.
pub fn render(particles: &[Particle], name: &str) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "TITLE     {name:<70}");
    for particle in particles {
        // Columns after the temperature factor (element/charge) are unused by the
        // CG model and left blank.
        let _ = writeln!(
            s,
            "ATOM  {:>5} {:>4} {:<4}A{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}{:>10}{:>2}{:>2} ",
            particle.atom_id % 100000,
            particle.bead,
            particle.residue,
            particle.molecule_id % 10000,
            particle.position[0] * ANGSTROM_PER_NM,
            particle.position[1] * ANGSTROM_PER_NM,
            particle.position[2] * ANGSTROM_PER_NM,
            0.0,
            0.0,
            "",
            "",
            "",
        );
    }
    let _ = writeln!(s, "TER");
    s
}

impl Membrane {
    /// Render this membrane as a `.pdb` coordinate file.
    pub fn to_pdb(&self) -> String {
        render(&self.particles, &self.name)
    }
}
