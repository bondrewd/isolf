//! GROMACS coordinate (`.gro`) rendering of a [`Membrane`].
//!
//! Particles are grouped by lipid species (sharing the canonical order used by
//! the `.top`) and renumbered sequentially. Residue and atom numbers wrap at
//! 100000, per the five-column `.gro` convention. Velocities are written as
//! zeros.

use std::fmt;

use crate::membrane::Membrane;

/// Wraps a [`Membrane`] for `.gro` rendering via [`fmt::Display`].
pub struct Gro<'a> {
    membrane: &'a Membrane,
}

impl<'a> Gro<'a> {
    /// Wrap a membrane for rendering.
    pub fn new(membrane: &'a Membrane) -> Self {
        Self { membrane }
    }
}

impl Membrane {
    /// Render this membrane as a GROMACS `.gro` coordinate file.
    pub fn to_gro(&self) -> String {
        Gro::new(self).to_string()
    }
}

impl fmt::Display for Gro<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let membrane = self.membrane;
        writeln!(f, "{}, t = {:>16.3} ", membrane.name, 0.0)?;
        writeln!(f, "{:>12} ", membrane.particles.len())?;

        // Particles are already in the canonical species-grouped order, so the
        // stored molecule/atom ids (wrapped at 100000 per the .gro convention)
        // are written directly.
        for particle in &membrane.particles {
            writeln!(
                f,
                "{:>5}{:>5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4} ",
                particle.molecule_id % 100000,
                particle.residue,
                particle.bead,
                particle.atom_id % 100000,
                particle.position[0],
                particle.position[1],
                particle.position[2],
                0.0,
                0.0,
                0.0,
            )?;
        }

        writeln!(
            f,
            "{:>15.4}{:>15.4}{:>15.4} ",
            membrane.box_size[0], membrane.box_size[1], membrane.box_size[2]
        )?;
        writeln!(f)
    }
}
