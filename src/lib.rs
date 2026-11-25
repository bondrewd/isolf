//! Build coarse-grained starting structures for GENESIS molecular dynamics:
//! iSoLF (Implicit Solvent Lipid Force Field) membranes, vesicles, and tubules,
//! AICG2+ / Clementi Gō coarse-grained proteins (standalone or embedded), and
//! 1D bilayer umbrella-sampling / REUS campaigns.
//!
//! Where to start, by area:
//!
//! - **Lipids & force fields.** [`force_field`] is the typed iSoLF force field,
//!   including [`force_field::lipid_protein`] for tunable `lipid_protein.itp`
//!   non-bonded tables; [`itp`] renders `isolf.itp`. [`membrane`] lays out a
//!   flat bilayer from per-leaflet [`composition`]s, [`vesicle`] wraps one into
//!   a sphere, and [`tubule`] wraps one into a periodic cylinder.
//! - **Proteins.** [`protein`] reads all-atom PDB/mmCIF files (or builds
//!   single-bead canonical amino acids) and constructs AICG2+ or Clementi Gō
//!   coarse-grained topologies; [`embed`] embeds one or more proteins into
//!   membranes, vesicles, or tubules with lipid carving and obstacle-aware
//!   relaxation.
//! - **Umbrella sampling.** [`umbrella`] places amino-acid or peptide/protein
//!   probes across bilayer `z` windows without carving lipids, extracts
//!   `RESTR_CVS` from GENESIS logs (with REUS `.rem` demuxing), and writes
//!   campaign manifests for MBAR and parameter-optimization workflows.
//! - **Packing.** [`pack`] is a Lloyd density-equalization / de-clash relaxation
//!   engine for lipid coordinates (lateral for a bilayer, over a shell for a
//!   vesicle or tubule), with pinned-site obstacle support.
//! - **Output.** Coordinates and topology in `.gro` ([`gro`]), `.top` ([`top`]),
//!   `.pdb` ([`pdb`]), `.psf` ([`psf`]), `.crd` ([`crd`]), and `.cif` ([`cif`]);
//!   GENESIS run-control and MBAR files ([`inp`]); VMD visualization scripts
//!   ([`vmd`]). Typed errors live in [`error`].
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
pub mod embed;
pub mod error;
pub mod force_field;
pub mod gro;
pub mod inp;
pub mod itp;
pub mod membrane;
pub mod pack;
pub mod pdb;
pub mod protein;
pub mod psf;
pub mod top;
pub mod tubule;
pub mod umbrella;
pub mod vesicle;
pub mod vmd;
pub mod voronoi;
