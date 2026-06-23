//! End-to-end checks of the membrane builder and its output writers.

use isolf::composition::Composition;
use isolf::force_field::ForceField;
use isolf::membrane::{BuildOptions, DEFAULT_NAME, Membrane, Sizing};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn pure(lipid: &str) -> Composition {
    Composition::from_weights([(lipid.to_string(), 1.0)]).unwrap()
}

fn small_membrane(lipid: &str, count: usize, seed: u64) -> Membrane {
    let force_field = ForceField::isolf();
    let leaflet = pure(lipid);
    let mut rng = StdRng::seed_from_u64(seed);
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
fn gro_has_header_count_atoms_and_box() {
    // 9 per leaflet ⇒ 18 DLPC lipids ⇒ 90 beads. (A multi-column grid so the
    // x/y box extent is non-zero.)
    let membrane = small_membrane("DLPC", 9, 1);
    let gro = membrane.to_gro();
    let lines: Vec<&str> = gro.lines().collect();

    assert!(lines[0].starts_with("CG membrane model, t ="));
    assert_eq!(lines[1].trim(), "90");

    // title + count, then the box and a trailing blank line, frame the atoms.
    let atom_lines = &lines[2..lines.len() - 2];
    assert_eq!(atom_lines.len(), 90);
    // Atom rows carry the (uppercase) residue and bead names.
    assert!(lines[2].contains("DLPC"));
    assert!(lines[2].contains("CHO"));

    // Box line: three whitespace-separated positive numbers.
    let box_line = lines[lines.len() - 2];
    let box_values: Vec<f64> = box_line
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(box_values.len(), 3);
    assert!(box_values.iter().all(|&v| v > 0.0));
}

#[test]
fn pdb_has_title_one_atom_per_bead_and_terminator() {
    let membrane = small_membrane("DLPC", 1, 1);
    let pdb = membrane.to_pdb();
    let lines: Vec<&str> = pdb.lines().collect();

    assert!(lines[0].starts_with("TITLE     CG membrane model"));
    assert_eq!(lines.iter().filter(|l| l.starts_with("ATOM")).count(), 10);
    assert_eq!(*lines.last().unwrap(), "TER");
    // Coordinates are in ångström (nm × 10): the head bead sits a few nm up,
    // i.e. tens of ångström.
    let first_atom = lines[1];
    assert!(first_atom.starts_with("ATOM      1  CHO DLPCA"));
}

#[test]
fn top_lists_molecules_and_chain_pairs() {
    let membrane = small_membrane("DLPC", 1, 1);
    let expected = "\
; Implicit Solvent Lipid Forcefield (iSoLF)
#include \"./isolf.itp\"

[ system ]
CG membrane model

[ molecules ]
DLPC 2

[ cg_ele_chain_pairs ]
ON 1 - 2 : 1 - 2
";
    assert_eq!(membrane.to_top(), expected);
}

#[test]
fn mixed_leaflet_counts_appear_in_top() {
    let force_field = ForceField::isolf();
    let leaflet = Composition::from_weights([("DOPC".into(), 1.0), ("DOPS".into(), 1.0)]).unwrap();
    let mut rng = StdRng::seed_from_u64(99);
    let membrane = Membrane::build(
        &force_field,
        DEFAULT_NAME,
        &leaflet,
        &leaflet,
        Sizing::Count {
            upper: 64,
            lower: 64,
        },
        &BuildOptions::default(),
        &mut rng,
    )
    .unwrap();

    let top = membrane.to_top();
    // 64 per leaflet, 50/50, doubled across leaflets.
    assert!(top.contains("DOPC 64"));
    assert!(top.contains("DOPS 64"));
    assert!(top.contains("ON 1 - 128 : 1 - 128"));
    // The .gro renders every bead.
    assert_eq!(
        membrane
            .to_gro()
            .lines()
            .filter(|l| l.contains("DOPC") || l.contains("DOPS"))
            .count(),
        membrane.particles.len()
    );
}

#[test]
fn placement_is_reproducible_for_a_fixed_seed() {
    let a = small_membrane("POPC", 36, 2024).to_gro();
    let b = small_membrane("POPC", 36, 2024).to_gro();
    assert_eq!(a, b);
}
