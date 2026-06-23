//! Build coarse-grained starting structures for GENESIS molecular dynamics:
//! iSoLF (Implicit Solvent Lipid Force Field) lipid membranes and vesicles.
//!
//! Where to start, by area:
//!
//! - **Lipids.** [`force_field`] is the typed iSoLF force field; [`itp`] renders
//!   it to a GENESIS/GROMACS `.itp`. [`membrane`] lays out a flat bilayer from
//!   per-leaflet [`composition`]s, and [`vesicle`] wraps one into a sphere.
//! - **Packing.** [`pack`] is a Lloyd density-equalization / de-clash relaxation
//!   engine for lipid coordinates (lateral for a bilayer, over a shell for a
//!   vesicle).
//! - **Output.** Coordinates and topology in `.gro` ([`gro`]), `.top` ([`top`]),
//!   `.pdb` ([`pdb`]), `.psf` ([`psf`]), `.crd` ([`crd`]), and `.cif` ([`cif`]);
//!   GENESIS run-control files ([`inp`]); VMD visualization scripts ([`vmd`]).
//!   Typed errors live in [`error`].
//!
//! ```
//! use isolf::composition::Composition;
//! use isolf::force_field::ForceField;
//! use isolf::membrane::{BuildOptions, Membrane, Sizing, DEFAULT_NAME};
//! use rand::SeedableRng;
//!
//! let force_field = ForceField::isolf();
//! let leaflet = Composition::from_weights([("POPC".to_string(), 1.0)]).unwrap();
//! let mut rng = rand::rngs::StdRng::seed_from_u64(0);
//! let membrane = Membrane::build(
//!     &force_field,
//!     DEFAULT_NAME,
//!     &leaflet,
//!     &leaflet,
//!     Sizing::Count { upper: 16, lower: 16 },
//!     &BuildOptions::default(),
//!     &mut rng,
//! )
//! .unwrap();
//! assert_eq!(membrane.total_lipids(), 32);
//! ```

pub mod cif;
pub mod composition;
pub mod crd;
pub mod error;
pub mod force_field;
pub mod gro;
pub mod inp;
pub mod itp;
pub mod membrane;
pub mod pack;
pub mod pdb;
pub mod psf;
pub mod top;
pub mod vesicle;
pub mod vmd;
pub mod voronoi;
