//! GROMACS coordinate (`.gro`) rendering of a [`Membrane`].
//!
//! Particles are grouped by lipid species (sharing the canonical order used by
//! the `.top`) and renumbered sequentially. Residue and atom numbers wrap at
//! 100000, per the five-column `.gro` convention. Velocities are written as
//! zeros.

use std::fmt::{self, Write as _};

use crate::error::BuildError;
use crate::force_field::ForceField;
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

pub(crate) fn write_gro_line(
    f: &mut impl fmt::Write,
    res_id: usize,
    residue: &str,
    bead: &str,
    atom_id: usize,
    position: [f64; 3],
) -> fmt::Result {
    writeln!(
        f,
        "{:>5}{:>5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}  0.0000  0.0000  0.0000 ",
        res_id % 100000,
        residue,
        bead,
        atom_id % 100000,
        position[0],
        position[1],
        position[2],
    )
}

fn write_gro_particle(f: &mut impl fmt::Write, p: &Particle) -> fmt::Result {
    write_gro_line(f, p.molecule_id, &p.residue, &p.bead, p.atom_id, p.position)
}

/// Render a particle list as a GROMACS `.gro` coordinate file using parallel
/// formatting for large systems.
pub fn render(particles: &[Particle], title: &str, box_size: [f64; 3]) -> String {
    let mut s = String::with_capacity(particles.len() * 70 + 128);
    let _ = writeln!(s, "{title}, t = {:>16.3} ", 0.0);
    let _ = writeln!(s, "{:>12} ", particles.len());
    write_particles(&mut s, particles, 70, |buf, p| {
        let _ = write_gro_particle(buf, p);
    });
    let _ = writeln!(
        s,
        "{:>15.4}{:>15.4}{:>15.4} ",
        box_size[0], box_size[1], box_size[2]
    );
    let _ = writeln!(s);
    s
}

impl Membrane {
    /// Render this membrane as a GROMACS `.gro` coordinate file.
    pub fn to_gro(&self) -> String {
        render(&self.particles, &self.name, self.box_size)
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

/// Parse an equilibrated iSoLF lipid bilayer `.gro` file into a [`Membrane`].
///
/// Bead charges and masses are resolved from `force_field`, and lipid molecule
/// counts are accumulated in order of appearance.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] on malformed `.gro` lines or unknown
/// lipid/bead names.
pub fn parse_bilayer(content: &str, force_field: &ForceField) -> Result<Membrane, BuildError> {
    let mut lines = content.lines();
    let raw_title = lines
        .next()
        .ok_or_else(|| BuildError::StructureParse("empty .gro file".to_string()))?;
    let name = raw_title
        .split(", t =")
        .next()
        .unwrap_or(raw_title)
        .trim()
        .to_string();
    let count_line = lines
        .next()
        .ok_or_else(|| BuildError::StructureParse("missing atom count line in .gro".to_string()))?;
    let n_atoms: usize = count_line.trim().parse().map_err(|_| {
        BuildError::StructureParse(format!("invalid .gro atom count: '{}'", count_line.trim()))
    })?;
    if n_atoms == 0 {
        return Err(BuildError::StructureParse(
            ".gro bilayer contains zero atoms".to_string(),
        ));
    }

    let mut particles = Vec::with_capacity(n_atoms);
    let mut cur_mol_id = 0usize;
    let mut prev_raw_res: Option<(usize, String)> = None;
    for idx in 0..n_atoms {
        let line = lines.next().ok_or_else(|| {
            BuildError::StructureParse(format!(".gro ended at atom {idx} of expected {n_atoms}"))
        })?;
        let (raw_res_id, res_name, bead_name, pos) = parse_gro_atom_line(line, idx + 1)?;
        let bt = force_field.bead_type(&bead_name).ok_or_else(|| {
            BuildError::StructureParse(format!(
                ".gro line {}: unknown iSoLF bead type '{bead_name}' (residue '{res_name}')",
                idx + 3
            ))
        })?;
        if force_field.lipid(&res_name).is_none() {
            return Err(BuildError::StructureParse(format!(
                ".gro line {}: unknown iSoLF lipid residue '{res_name}'",
                idx + 3
            )));
        }
        let is_new_mol = match &prev_raw_res {
            Some((prev_id, prev_res)) => *prev_id != raw_res_id || *prev_res != res_name,
            None => true,
        };
        if is_new_mol {
            cur_mol_id += 1;
            prev_raw_res = Some((raw_res_id, res_name.clone()));
        }
        particles.push(Particle {
            residue: res_name,
            bead: bead_name,
            molecule_id: cur_mol_id,
            atom_id: idx + 1,
            position: pos,
            charge: bt.charge,
            mass: bt.mass,
        });
    }

    let box_line = lines
        .find(|l| !l.trim().is_empty())
        .ok_or_else(|| BuildError::StructureParse("missing box vector line in .gro".to_string()))?;
    let box_vals: Vec<f64> = box_line
        .split_whitespace()
        .take(3)
        .map(|tok| tok.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| {
            BuildError::StructureParse(format!("invalid .gro box line: '{}'", box_line.trim()))
        })?;
    if box_vals.len() < 3 || box_vals.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(BuildError::StructureParse(format!(
            "invalid .gro box dimensions: '{}'",
            box_line.trim()
        )));
    }
    let box_size = [box_vals[0], box_vals[1], box_vals[2]];

    // Preserve contiguous molecule runs in exact coordinate order so `[ molecules ]`
    // in the generated topology matches interleaved `.gro` bilayers atom-for-atom.
    let mut lipid_counts: Vec<(String, usize)> = Vec::new();
    for mol in particles.chunk_by(|a, b| a.molecule_id == b.molecule_id) {
        let res = &mol[0].residue;
        if let Some((last_res, count)) = lipid_counts.last_mut()
            && last_res == res
        {
            *count += 1;
        } else {
            lipid_counts.push((res.clone(), 1));
        }
    }

    Ok(Membrane {
        name,
        box_size,
        lipid_counts,
        particles,
    })
}

fn parse_gro_atom_line(
    line: &str,
    atom_num: usize,
) -> Result<(usize, String, String, [f64; 3]), BuildError> {
    if line.len() < 44 {
        return Err(BuildError::StructureParse(format!(
            ".gro atom {atom_num} line is too short ({} chars)",
            line.len()
        )));
    }
    let res_id: usize = line[0..5].trim().parse().map_err(|_| {
        BuildError::StructureParse(format!(
            ".gro atom {atom_num}: invalid residue number '{}'",
            line[0..5].trim()
        ))
    })?;
    let res_name = line[5..10].trim().to_ascii_uppercase();
    let bead_name = line[10..15].trim().to_ascii_uppercase();

    // Standard 8-char columns (20..28, 28..36, 36..44), with whitespace fallback
    // for wider coordinate columns.
    let coords = if let (Ok(x), Ok(y), Ok(z)) = (
        line[20..28].trim().parse::<f64>(),
        line[28..36].trim().parse::<f64>(),
        line[36..44].trim().parse::<f64>(),
    ) {
        [x, y, z]
    } else {
        let nums: Vec<f64> = line[20..]
            .split_whitespace()
            .take(3)
            .map(|t| t.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| {
                BuildError::StructureParse(format!(
                    ".gro atom {atom_num}: invalid coordinates in '{}'",
                    line.trim()
                ))
            })?;
        if nums.len() < 3 {
            return Err(BuildError::StructureParse(format!(
                ".gro atom {atom_num}: expected 3 coordinates"
            )));
        }
        [nums[0], nums[1], nums[2]]
    };
    if coords.iter().any(|c| !c.is_finite()) {
        return Err(BuildError::StructureParse(format!(
            ".gro atom {atom_num}: non-finite coordinates"
        )));
    }
    Ok((res_id, res_name, bead_name, coords))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composition::Composition;
    use crate::membrane::{BuildOptions, Sizing};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn parse_bilayer_round_trips_generated_membrane_and_preserves_interleaved_runs() {
        let ff = ForceField::isolf();
        let comp = Composition::parse("POPC=1,POPS=1").unwrap();
        let mut rng = StdRng::seed_from_u64(3);
        let mem = Membrane::build(
            &ff,
            "mixed",
            &comp,
            &comp,
            Sizing::Count { upper: 8, lower: 8 },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap();

        let gro_text = mem.to_gro();
        let parsed = parse_bilayer(&gro_text, &ff).unwrap();
        assert_eq!(parsed.particles.len(), mem.particles.len());
        assert_eq!(parsed.lipid_counts, mem.lipid_counts);
        for (a, b) in parsed.box_size.iter().zip(&mem.box_size) {
            assert!((a - b).abs() < 1e-3);
        }

        // Interleaved species in an external .gro must preserve contiguous runs so
        // `[ molecules ]` stays aligned atom-for-atom with coordinates.
        let popc = ff.lipid("POPC").unwrap();
        let pops = ff.lipid("POPS").unwrap();
        let mut parts = Vec::new();
        let mut atom_id = 0usize;
        for (mol_id, (res, spec)) in [
            (1, ("POPC", popc)),
            (2, ("POPS", pops)),
            (3, ("POPC", popc)),
        ] {
            for bead in &spec.beads {
                atom_id += 1;
                let bt = ff.bead_type(bead).unwrap();
                parts.push(Particle {
                    residue: res.to_string(),
                    bead: bead.clone(),
                    molecule_id: mol_id,
                    atom_id,
                    position: [0.1 * atom_id as f64, 0.2, 0.3],
                    charge: bt.charge,
                    mass: bt.mass,
                });
            }
        }
        let interleaved = render(&parts, "interleaved", [8.0, 8.0, 20.0]);
        let parsed_interleaved = parse_bilayer(&interleaved, &ff).unwrap();
        assert_eq!(
            parsed_interleaved.lipid_counts,
            vec![
                ("POPC".to_string(), 1),
                ("POPS".to_string(), 1),
                ("POPC".to_string(), 1),
            ]
        );
    }
}
