//! PDB coordinate (`.pdb`) rendering of a particle list.
//!
//! Particles are written in order with their molecule and atom numbers.
//! Coordinates are converted from nm to ångström (×10), as the PDB format
//! expects. A [`Membrane`] renders through the [`render`] core.

use std::fmt::Write;

use crate::membrane::{Membrane, Particle, write_particles};

/// nm → ångström conversion for PDB coordinates.
const ANGSTROM_PER_NM: f64 = 10.0;

/// Append one PDB `ATOM` record (coordinates converted from nm to Å).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_atom_line(
    buf: &mut String,
    atom_id: usize,
    bead: &str,
    residue: &str,
    chain: char,
    res_id: usize,
    position: [f64; 3],
    occ: &str,
) {
    // Columns after the temperature factor (element/charge) are unused by the
    // CG model and left blank.
    let _ = writeln!(
        buf,
        "ATOM  {:>5} {:>4} {:<4}{}{:>4}    {:>8.3}{:>8.3}{:>8.3}  {occ}  0.00               ",
        atom_id % 100000,
        bead,
        residue,
        chain,
        res_id % 10000,
        position[0] * ANGSTROM_PER_NM,
        position[1] * ANGSTROM_PER_NM,
        position[2] * ANGSTROM_PER_NM,
    );
}

/// Render a particle list as a `.pdb` coordinate file.
pub fn render(particles: &[Particle], name: &str) -> String {
    let mut s = String::with_capacity(particles.len() * 82 + 96);
    let _ = writeln!(s, "TITLE     {name:<70}");
    write_particles(&mut s, particles, 82, |buf, p| {
        write_atom_line(
            buf,
            p.atom_id,
            &p.bead,
            &p.residue,
            'A',
            p.molecule_id,
            p.position,
            "0.00",
        );
    });
    let _ = writeln!(s, "TER");
    s
}

impl Membrane {
    /// Render this membrane as a `.pdb` coordinate file.
    pub fn to_pdb(&self) -> String {
        render(&self.particles, &self.name)
    }
}
