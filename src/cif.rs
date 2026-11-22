//! mmCIF (`.cif`) rendering of a particle list.
//!
//! Emits a minimal `atom_site` loop with cartesian coordinates in ångström
//! (nm × 10). Coarse-grained beads have no real chemical element, so every atom
//! is given a placeholder `type_symbol` of `C` for the benefit of structure
//! viewers. A [`Membrane`] renders through the [`render`] core.

use std::fmt::Write;

use crate::membrane::{Membrane, Particle, write_particles};

/// nm → ångström conversion for cartesian coordinates.
const ANGSTROM_PER_NM: f64 = 10.0;

/// Render a particle list as an mmCIF `.cif` file.
pub fn render(particles: &[Particle], name: &str) -> String {
    let mut s = String::with_capacity(particles.len() * 64 + 320);
    let _ = writeln!(s, "data_{}", data_block_name(name));
    let _ = writeln!(s, "#");
    let _ = writeln!(s, "loop_");
    for column in [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
    ] {
        let _ = writeln!(s, "{column}");
    }

    write_particles(&mut s, particles, 64, |buf, particle| {
        let [x, y, z] = particle.position.map(|c| c * ANGSTROM_PER_NM);
        let _ = writeln!(
            buf,
            "ATOM {} C {} {} A {} {:.3} {:.3} {:.3} 1.00 0.00",
            particle.atom_id, particle.bead, particle.residue, particle.molecule_id, x, y, z,
        );
    });
    let _ = writeln!(s, "#");
    s
}

impl Membrane {
    /// Render this membrane as an mmCIF `.cif` file.
    pub fn to_cif(&self) -> String {
        render(&self.particles, &self.name)
    }
}

/// Turn the system name into a whitespace-free mmCIF data-block identifier.
fn data_block_name(name: &str) -> String {
    let sanitized: String = name.split_whitespace().collect::<Vec<_>>().join("_");
    if sanitized.is_empty() {
        "membrane".to_string()
    } else {
        sanitized
    }
}

#[cfg(test)]
mod tests {
    use crate::composition::Composition;
    use crate::force_field::ForceField;
    use crate::membrane::{BuildOptions, DEFAULT_NAME, Membrane, Sizing};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn membrane() -> Membrane {
        let force_field = ForceField::isolf();
        let leaflet = Composition::parse("DLPC=1").unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        Membrane::build(
            &force_field,
            DEFAULT_NAME,
            &leaflet,
            &leaflet,
            Sizing::Count { upper: 4, lower: 4 },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap()
    }

    #[test]
    fn has_atom_site_loop_and_all_atoms() {
        let membrane = membrane();
        let cif = membrane.to_cif();

        assert!(cif.starts_with("data_CG_membrane_model\n"));
        assert!(cif.contains("loop_"));
        assert!(cif.contains("_atom_site.Cartn_x"));
        // 8 DLPC × 5 beads = 40 ATOM records.
        assert_eq!(cif.lines().filter(|l| l.starts_with("ATOM ")).count(), 40);
        assert!(cif.trim_end().ends_with('#'));
    }
}
