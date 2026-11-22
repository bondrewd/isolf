//! GROMACS coordinate (`.gro`) rendering of a [`Membrane`].
//!
//! Particles are grouped by lipid species (sharing the canonical order used by
//! the `.top`) and renumbered sequentially. Residue and atom numbers wrap at
//! 100000, per the five-column `.gro` convention. Velocities are written as
//! zeros.

use std::fmt::{self, Write as _};

use crate::membrane::{Membrane, Particle, write_particles};

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

fn write_gro_particle(f: &mut impl fmt::Write, p: &Particle) -> fmt::Result {
    writeln!(
        f,
        "{:>5}{:>5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}  0.0000  0.0000  0.0000 ",
        p.molecule_id % 100000,
        p.residue,
        p.bead,
        p.atom_id % 100000,
        p.position[0],
        p.position[1],
        p.position[2],
    )
}

impl Membrane {
    /// Render this membrane as a GROMACS `.gro` coordinate file.
    pub fn to_gro(&self) -> String {
        let mut s = String::with_capacity(self.particles.len() * 70 + 128);
        let _ = writeln!(s, "{}, t = {:>16.3} ", self.name, 0.0);
        let _ = writeln!(s, "{:>12} ", self.particles.len());
        write_particles(&mut s, &self.particles, 70, |buf, p| {
            let _ = write_gro_particle(buf, p);
        });
        let _ = writeln!(
            s,
            "{:>15.4}{:>15.4}{:>15.4} ",
            self.box_size[0], self.box_size[1], self.box_size[2]
        );
        let _ = writeln!(s);
        s
    }
}

impl fmt::Display for Gro<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let membrane = self.membrane;
        writeln!(f, "{}, t = {:>16.3} ", membrane.name, 0.0)?;
        writeln!(f, "{:>12} ", membrane.particles.len())?;
        for particle in &membrane.particles {
            write_gro_particle(f, particle)?;
        }
        writeln!(
            f,
            "{:>15.4}{:>15.4}{:>15.4} ",
            membrane.box_size[0], membrane.box_size[1], membrane.box_size[2]
        )?;
        writeln!(f)
    }
}
