//! Command-line builder for coarse-grained iSoLF membranes, vesicles, tubules,
//! and AICG2+ proteins for GENESIS.
//!
//! ```sh
//! isolf membrane --lipids "POPC=3,POPS=1" --count 512 --inp --vmd -o ./sim
//! isolf membrane --lipids POPC --box 20 -p channel.pdb --inp --vmd -o ./prot_mem
//! isolf vesicle  --lipids POPC --radius 20 --top -o ./ves
//! isolf tubule   --lipids "DPPC=1,POPC=3" --phase raft --radius 10 --length 25 -o ./tube
//! isolf protein  1aki.pdb --inp --vmd -o ./prot
//! isolf lipids
//! ```
//!
//! `<name>.gro` and `<name>.log` are always written. `--top` adds the topology
//! (`<name>.top` plus `isolf.itp`, and `aicg2p.itp` + `<prot>.itp` when proteins
//! are present); `--inp` (or any `--min`/`--equil`/`--prod` override) adds
//! GENESIS control files and implies `--top`; `--vmd` (or any `--vmd-*` flag)
//! adds VMD scripts and implies `--psf`; `--gif` (or any `--gif-*` flag) adds a
//! relaxation animation; `.pdb`, `.psf`, `.crd`, and `.cif` are added on request.

use std::error::Error;
use std::f64::consts::TAU;
use std::fmt::Write as _;
use std::fs;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};

use clap::builder::styling::{AnsiColor, Styles};
use clap::{ColorChoice, CommandFactory, FromArgMatches, Parser, Subcommand};
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;

use isolf::composition::Composition;
use isolf::embed::{
    EmbeddedSystem, ProteinPlacement, ProteinSpec, ZPlacement, embed_membrane,
    embed_membrane_recorded, embed_tubule, embed_tubule_recorded, embed_vesicle,
    embed_vesicle_recorded,
};
use isolf::error::BuildError;
use isolf::force_field::ForceField;
use isolf::membrane::{
    BuildOptions, DEFAULT_NAME, DEFAULT_PADDING, DEFAULT_TEMPERATURE, Phase, RelaxFrame, Sizing,
};
use isolf::protein::Protein;
use isolf::protein::aicg2p::{
    AicgScale, ProteinModel, ProteinOptions, SafeDihedral, parse_charge_table, parse_selection,
};
use isolf::protein::params::render_aicg2p_itp;
use isolf::tubule::TubuleRadius;
use isolf::vesicle::VesicleRadius;
use isolf::vmd::{Layout, VmdOptions, VmdSource};

mod anim;
mod report;
mod runlog;
mod update;

/// Where the finished system is centred (`--center`). The periodic box is the
/// same either way; only the coordinates move.
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum CenterArg {
    /// Centred on the coordinate origin (0, 0, 0) (default).
    #[default]
    Origin,
    /// Centred in the periodic box [0, box].
    Box,
}

/// Scope for `--vmd-color-by`: how lipid beads are coloured. Maps to
/// [`isolf::vmd::LipidColoring`].
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum LipidColorScope {
    /// Colour the head bead by lipid species; the rest by bead role.
    Head,
    /// Colour every bead by its role (head/phosphate/glycerol/tail).
    Bead,
    /// Colour every bead by lipid species.
    Whole,
}

/// Scope for `--gif-mode`: how `--gif` renders each frame. Maps to
/// [`anim::GifMode`].
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum GifModeArg {
    /// One disc per lipid; clear for small to medium systems.
    Point,
    /// Local-density heatmap; better for large, dense systems.
    Density,
}

/// Initial lateral lipid distribution (`--phase`). Maps to [`Phase`].
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum PhaseArg {
    /// Uniformly mixed lipids (default).
    #[default]
    Random,
    /// Parallel stripes (membrane), latitude bands (vesicle), or axial rings (tubule).
    Stripe,
    /// Circular raft / polar cap for the first domain in a mixed matrix.
    Raft,
}

impl From<PhaseArg> for Phase {
    fn from(phase: PhaseArg) -> Self {
        match phase {
            PhaseArg::Random => Phase::Random,
            PhaseArg::Stripe => Phase::Stripe,
            PhaseArg::Raft => Phase::Raft,
        }
    }
}

/// Where a standalone protein is centred (`isolf protein --center`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum ProteinCenterArg {
    /// Keep the input PDB/mmCIF coordinates unchanged (default).
    #[default]
    Keep,
    /// Translate centroid to the coordinate origin (0, 0, 0).
    Origin,
    /// Translate centroid to the centre of the periodic box [0, box].
    Box,
}

/// Protein coarse-grained model (`--model`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum ProteinModelArg {
    /// AICG2+ structure-based model with flexible local potentials (default).
    #[default]
    #[value(name = "aicg2+", alias = "aicg2p", alias = "aicg2")]
    Aicg2Plus,
    /// Clementi Gō-like Cα model.
    #[value(name = "clementi", alias = "ccgo")]
    ClementiGo,
}

impl From<ProteinModelArg> for ProteinModel {
    fn from(m: ProteinModelArg) -> Self {
        match m {
            ProteinModelArg::Aicg2Plus => ProteinModel::Aicg2Plus,
            ProteinModelArg::ClementiGo => ProteinModel::ClementiGo,
        }
    }
}

/// AICG2+ energy scaling scheme (`--scale`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum AicgScaleArg {
    /// Universal AICG2+ scaling factors (1.11 / 0.87 / 0.32, default).
    #[default]
    General,
    /// Normalize mean local and contact energies to reference averages (1.72 / 1.23 / 0.55).
    Average,
}

impl From<AicgScaleArg> for AicgScale {
    fn from(s: AicgScaleArg) -> Self {
        match s {
            AicgScaleArg::General => AicgScale::General,
            AicgScaleArg::Average => AicgScale::Average,
        }
    }
}

/// Safe dihedral damping mode (`--safe-dihedral`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum SafeDihedralArg {
    /// Unmodified dihedral potentials (0).
    #[value(name = "0", alias = "none")]
    None,
    /// cos²(kθ)-damped safe dihedrals (1, default).
    #[default]
    #[value(name = "1", alias = "cos2")]
    Cos2,
    /// sin³(kθ)-damped safe dihedrals (3).
    #[value(name = "3", alias = "sin3")]
    Sin3,
}

impl From<SafeDihedralArg> for SafeDihedral {
    fn from(s: SafeDihedralArg) -> Self {
        match s {
            SafeDihedralArg::None => SafeDihedral::None,
            SafeDihedralArg::Cos2 => SafeDihedral::Cos2,
            SafeDihedralArg::Sin3 => SafeDihedral::Sin3,
        }
    }
}

/// Build coarse-grained iSoLF membranes, vesicles, tubules, and AICG2+ proteins for GENESIS.
#[derive(Parser, Debug)]
#[command(
    name = "isolf",
    version,
    about,
    disable_help_subcommand = true,
    arg_required_else_help = true
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build a flat lipid bilayer.
    Membrane(MembraneArgs),
    /// Build a spherical lipid vesicle.
    Vesicle(VesicleArgs),
    /// Build a cylindrical lipid tubule (periodic along z).
    Tubule(TubuleArgs),
    /// Coarse-grain a protein PDB/mmCIF with AICG2+ or Clementi Gō for GENESIS.
    Protein(ProteinArgs),
    /// List all supported iSoLF lipids.
    Lipids(LipidsArgs),
    /// Update isolf in place to the latest release.
    Update,
    /// Uninstall isolf, removing the binary from the system.
    Uninstall,
}

/// Arguments for `isolf membrane`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Build a flat lipid bilayer",
    long_about = "\
Build a flat coarse-grained lipid bilayer for GENESIS.

Size the bilayer with --count (lipid count per leaflet), --box (in-plane nm),
or both (pack an exact count into a fixed box). Use --lipids for a symmetric
bilayer, or --upper / --lower for asymmetric leaflets."
)]
struct MembraneArgs {
    // ---- Lipid composition ----
    /// Lipid composition for both leaflets (symmetric bilayer) as NAME=WEIGHT pairs.
    ///
    /// Examples:
    ///   -L POPC
    ///   -L "POPC=3,POPS=1"
    ///   -L "POPC=1/2,DOPS=1/2"
    ///   -L "POPC=70%,POPS=30%"
    #[arg(
        short = 'L',
        long,
        value_name = "SPEC",
        conflicts_with_all = ["upper", "lower"],
        verbatim_doc_comment,
        help_heading = "Lipid composition"
    )]
    lipids: Vec<String>,

    /// Upper leaflet composition as NAME=WEIGHT pairs (mirrored to lower if --lower is omitted).
    #[arg(
        short = 'u',
        long,
        value_name = "SPEC",
        help_heading = "Lipid composition"
    )]
    upper: Vec<String>,

    /// Lower leaflet composition as NAME=WEIGHT pairs (mirrored to upper if --upper is omitted).
    #[arg(
        short = 'l',
        long,
        value_name = "SPEC",
        help_heading = "Lipid composition"
    )]
    lower: Vec<String>,

    // ---- Geometry ----
    /// Lipids per leaflet: <N> for both, or "up=<N>,lo=<N>".
    ///
    /// Examples:
    ///   -n 256                 # 256 lipids per leaflet in a square box
    ///   -n "up=300,lo=200"     # different counts per leaflet
    ///   -n 256 -b 20           # pack 256 lipids per leaflet into a 20 nm box
    #[arg(
        short = 'n',
        long,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    count: Vec<String>,

    /// In-plane box size in nm: <L> for a square, or "x=<L>,y=<L>".
    ///
    /// Examples:
    ///   -b 20                  # 20 x 20 nm square membrane
    ///   -b "x=30,y=15"         # 30 x 15 nm rectangular membrane
    #[arg(
        short = 'b',
        long = "box",
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    box_size: Vec<String>,

    /// System centre: origin (0,0,0) or box (Lx/2,Ly/2,Lz/2).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = CenterArg::Origin,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    center: CenterArg,

    /// Vacuum padding in nm added on each side along z.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    padding: f64,

    /// Keep unequal leaflet counts as given instead of growing the lighter leaflet to match area.
    #[arg(long, help_heading = "Geometry")]
    no_balance: bool,

    #[command(flatten)]
    common: CommonBuildArgs,
}

/// Arguments for `isolf vesicle`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Build a spherical lipid vesicle",
    long_about = "\
Build a spherical coarse-grained lipid vesicle for GENESIS.

Set the size with --radius: a bare number or ro=<R> fixes the outer radius,
ri=<R> fixes the inner (lumen) radius, and ri=<R>,ro=<R> fixes both. Use
--lipids for a symmetric vesicle, or --outer / --inner for asymmetric shells."
)]
struct VesicleArgs {
    // ---- Lipid composition ----
    /// Lipid composition for both shells (symmetric vesicle) as NAME=WEIGHT pairs.
    #[arg(
        short = 'L',
        long,
        value_name = "SPEC",
        conflicts_with_all = ["outer", "inner"],
        help_heading = "Lipid composition"
    )]
    lipids: Vec<String>,

    /// Outer shell composition as NAME=WEIGHT pairs (copied to inner if --inner is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    outer: Vec<String>,

    /// Inner (lumen) shell composition as NAME=WEIGHT pairs (copied to outer if --outer is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    inner: Vec<String>,

    // ---- Geometry ----
    /// Vesicle radius in nm: <R> or "ro=<R>" (outer), "ri=<R>" (inner), or both.
    ///
    /// Examples:
    ///   -r 20                  # outer radius 20 nm, inner derived automatically
    ///   -r "ri=15"             # inner radius 15 nm, outer derived automatically
    ///   -r "ri=14,ro=20"       # both radii fixed
    #[arg(
        short = 'r',
        long,
        required = true,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    radius: Vec<String>,

    /// System centre: origin (0,0,0) or box (Lx/2,Ly/2,Lz/2).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = CenterArg::Origin,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    center: CenterArg,

    /// Vacuum padding in nm added on every side of the cubic box.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    padding: f64,

    #[command(flatten)]
    common: CommonBuildArgs,
}

/// Arguments for `isolf tubule`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Build a cylindrical lipid tubule (periodic along z)",
    long_about = "\
Build a cylindrical coarse-grained lipid tubule periodic along z for GENESIS.

Set the cross-section with --radius (a bare number or ro=<R> fixes the outer
radius, ri=<R> fixes the inner lumen radius, ri=<R>,ro=<R> fixes both) and the
periodic axial length with --length. Use --lipids for a symmetric tubule, or
--outer / --inner for asymmetric shells."
)]
struct TubuleArgs {
    // ---- Lipid composition ----
    /// Lipid composition for both shells (symmetric tubule) as NAME=WEIGHT pairs.
    #[arg(
        short = 'L',
        long,
        value_name = "SPEC",
        conflicts_with_all = ["outer", "inner"],
        help_heading = "Lipid composition"
    )]
    lipids: Vec<String>,

    /// Outer shell composition as NAME=WEIGHT pairs (copied to inner if --inner is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    outer: Vec<String>,

    /// Inner (lumen) shell composition as NAME=WEIGHT pairs (copied to outer if --outer is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    inner: Vec<String>,

    // ---- Geometry ----
    /// Tubule radius in nm: <R> or "ro=<R>" (outer), "ri=<R>" (inner), or both.
    ///
    /// Examples:
    ///   -r 10                  # outer radius 10 nm, inner derived automatically
    ///   -r "ri=6"              # inner radius 6 nm, outer derived automatically
    ///   -r "ri=6,ro=10.5"      # both radii fixed
    #[arg(
        short = 'r',
        long,
        required = true,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    radius: Vec<String>,

    /// Periodic axial length of the tubule along z in nm.
    #[arg(
        short = 'l',
        long,
        required = true,
        value_name = "NM",
        help_heading = "Geometry"
    )]
    length: f64,

    /// System centre: origin (0,0,0) or box (Lx/2,Ly/2,Lz/2).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = CenterArg::Origin,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    center: CenterArg,

    /// Vacuum padding in nm added on each side along x and y.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    padding: f64,

    #[command(flatten)]
    common: CommonBuildArgs,
}

/// Options shared between `membrane`, `vesicle`, and `tubule`.
#[derive(clap::Args, Debug)]
struct CommonBuildArgs {
    // ---- Lipid composition (phase separation) ----
    /// Lateral domain layout: random (mixed), stripe (bands/rings), or raft (circular domain).
    ///
    /// Separate explicit multi-lipid domains with ';', e.g.
    /// -L "DPPC=1,DSPC=1;POPC=2,POPS=1" --phase raft
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = PhaseArg::Random,
        value_name = "MODE",
        verbatim_doc_comment,
        help_heading = "Lipid composition"
    )]
    phase: PhaseArg,

    /// Offset phase domains on opposite leaflets/shells (anti-registered domains).
    #[arg(long, help_heading = "Lipid composition")]
    anti_register: bool,

    // ---- Proteins ----
    /// Embed a PDB/mmCIF protein (repeatable) as PATH[,KEY=VAL,...].
    ///
    /// Preserves input orientation by default (OPM midplane z = 0).
    /// Keys: z=<nm|center>, x=<nm>, y=<nm>, tilt=<deg>, rot=<deg>, copies=<N>,
    ///       clearance=<nm>, model=<aicg2+|clementi>, scale=<general|average>,
    ///       safe-dihedral=<0|1|3>, idr=<ranges>, idr-local=<ranges>,
    ///       idr-nonlocal=<ranges>, charges=<file>, name=<molname>
    /// Example: -p channel.pdb,z=center,copies=2,clearance=0.2
    #[arg(
        short = 'p',
        long = "protein",
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Proteins"
    )]
    proteins: Vec<String>,

    // ---- Output files ----
    /// Output directory.
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output files"
    )]
    out: PathBuf,

    /// Base name for output files [default: membrane / vesicle / tubule].
    #[arg(long, value_name = "NAME", help_heading = "Output files")]
    name: Option<String>,

    /// Write topology (.top) and force-field (isolf.itp) files.
    #[arg(long, help_heading = "Output files")]
    top: bool,

    /// Write GENESIS .inp control files for min, equil, and prod (implies --top).
    #[arg(long, help_heading = "Output files")]
    inp: bool,

    /// Write VMD .vmd visualization scripts (implies --psf).
    #[arg(long, help_heading = "Output files")]
    vmd: bool,

    /// Write a .gif animating the leaflet relaxation.
    #[arg(long, help_heading = "Output files")]
    gif: bool,

    /// Write a .pdb coordinate file.
    #[arg(long, help_heading = "Output files")]
    pdb: bool,

    /// Write a CHARMM .psf structure file.
    #[arg(long, help_heading = "Output files")]
    psf: bool,

    /// Write a CHARMM .crd coordinate file.
    #[arg(long, help_heading = "Output files")]
    crd: bool,

    /// Write an mmCIF .cif coordinate file.
    #[arg(long, help_heading = "Output files")]
    cif: bool,

    // ---- GENESIS control files (--inp) ----
    /// Simulation temperature in K for GENESIS .inp files [default: 323.15].
    #[arg(
        short = 't',
        long,
        value_name = "K",
        help_heading = "GENESIS control files (--inp)"
    )]
    temperature: Option<f64>,

    /// Minimization (min.inp) overrides as KEY=VAL or <steps> (implies --inp).
    ///
    /// Keys: steps (1000), eneout (10), crdout (1000), rstout (1000), nbupdate (10).
    /// Example: --min "steps=2000,eneout=20"
    #[arg(
        long,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "GENESIS control files (--inp)"
    )]
    min: Vec<String>,

    /// Equilibration (npt.inp / nvt.inp) overrides as KEY=VAL or <steps> (implies --inp).
    ///
    /// Keys: steps (20000), dt (0.01), temp (--temperature), eneout (100),
    ///       crdout (100), rstout (1000), stoptr (10), nbupdate (10), seed (random).
    /// Example: --equil "steps=50000,temp=310"
    #[arg(
        long,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "GENESIS control files (--inp)"
    )]
    equil: Vec<String>,

    /// Production (pro.inp) overrides as KEY=VAL or <steps> (implies --inp).
    ///
    /// Keys: steps (500000), dt (0.02), temp (--temperature), eneout (500),
    ///       crdout (500), rstout (1000), stoptr (10), nbupdate (10), seed (random).
    /// Example: --prod "steps=1000000,crdout=1000"
    #[arg(
        long,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "GENESIS control files (--inp)"
    )]
    prod: Vec<String>,

    // ---- VMD visualization (--vmd) ----
    /// Use a smooth QuickSurf surface instead of VDW spheres (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    vmd_surface: bool,

    /// Slice the system with a clipping plane to reveal the cross-section (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    vmd_cutaway: bool,

    /// Append a Tachyon ray-tracing render command to each script (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    vmd_render: bool,

    /// VMD background colour name [default: white] (implies --vmd).
    #[arg(long, value_name = "COLOR", help_heading = "VMD visualization (--vmd)")]
    vmd_background: Option<String>,

    /// VMD representation material, e.g. AOChalky or Goodsell [default: AOChalky] (implies --vmd).
    #[arg(long, value_name = "NAME", help_heading = "VMD visualization (--vmd)")]
    vmd_material: Option<String>,

    /// Trajectory smoothing window in frames [default: 5] (implies --vmd).
    #[arg(long, value_name = "N", help_heading = "VMD visualization (--vmd)")]
    vmd_smooth: Option<usize>,

    /// Colouring scheme [default: head] (implies --vmd).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        value_name = "SCOPE",
        help_heading = "VMD visualization (--vmd)"
    )]
    vmd_color_by: Option<LipidColorScope>,

    // ---- Relaxation animation (--gif) ----
    /// Render mode: point scatter or density heatmap [default: point] (implies --gif).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        value_name = "MODE",
        help_heading = "Relaxation animation (--gif)"
    )]
    gif_mode: Option<GifModeArg>,

    /// Resolution multiplier, e.g. 2 doubles it [default: 1] (implies --gif).
    #[arg(
        long,
        value_name = "FACTOR",
        help_heading = "Relaxation animation (--gif)"
    )]
    gif_scale: Option<f64>,

    /// Playback speed in frames per second [default: 14] (implies --gif).
    #[arg(
        long,
        value_name = "FPS",
        help_heading = "Relaxation animation (--gif)"
    )]
    gif_fps: Option<u16>,

    // ---- General ----
    /// Seed value for the RNG (random if unset).
    #[arg(long, value_name = "SEED", help_heading = "General")]
    seed: Option<u64>,

    /// Print only the final summary and any warnings.
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        help_heading = "General"
    )]
    quiet: bool,

    /// Print full file paths in the Write section.
    #[arg(short = 'v', long, help_heading = "General")]
    verbose: bool,

    /// Disable coloured terminal output.
    #[arg(long, help_heading = "General")]
    no_color: bool,

    /// Use ASCII status markers instead of Unicode glyphs.
    #[arg(long, help_heading = "General")]
    ascii: bool,
}

impl CommonBuildArgs {
    fn want_inp(&self) -> bool {
        self.inp || !self.min.is_empty() || !self.equil.is_empty() || !self.prod.is_empty()
    }

    fn want_vmd(&self) -> bool {
        self.vmd
            || self.vmd_surface
            || self.vmd_cutaway
            || self.vmd_render
            || self.vmd_background.is_some()
            || self.vmd_material.is_some()
            || self.vmd_smooth.is_some()
            || self.vmd_color_by.is_some()
    }

    fn want_gif(&self) -> bool {
        self.gif || self.gif_mode.is_some() || self.gif_scale.is_some() || self.gif_fps.is_some()
    }
}

/// Arguments for `isolf protein`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Coarse-grain a protein PDB/mmCIF with AICG2+ or Clementi Gō for GENESIS",
    long_about = "\
Coarse-grain an all-atom protein PDB or mmCIF structure using the AICG2+
(default) or Clementi Gō Cα model for GENESIS RESIDCG.

Always writes <name>.gro, <name>.itp, aicg2p.itp, <name>.top, and <name>.log.
Add --inp for GENESIS min/nvt/pro control files, --vmd for VMD visualization
scripts (implies --psf), or --pdb / --psf / --crd / --cif as needed."
)]
struct ProteinArgs {
    /// Input all-atom protein structure (.pdb or .cif/.mmcif).
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    // ---- Protein model ----
    /// Coarse-grained protein force field: aicg2+ (default) or clementi.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = ProteinModelArg::Aicg2Plus,
        value_name = "MODEL",
        help_heading = "Protein model"
    )]
    model: ProteinModelArg,

    /// AICG2+ energy scaling scheme: general (default) or average.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = AicgScaleArg::General,
        value_name = "SCHEME",
        help_heading = "Protein model"
    )]
    scale: AicgScaleArg,

    /// Safe dihedral potential mode: 0 (none), 1 (cos2, default), or 3 (sin3).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = SafeDihedralArg::Cos2,
        value_name = "MODE",
        help_heading = "Protein model"
    )]
    safe_dihedral: SafeDihedralArg,

    /// Disordered residue ranges excluded from both local and nonlocal Gō terms (e.g. "1-25,80-95").
    #[arg(long, value_name = "RANGES", help_heading = "Protein model")]
    idr: Option<String>,

    /// Disordered residue ranges excluded from AICG2+ local Gaussian 1-3 and 1-4 terms.
    #[arg(long, value_name = "RANGES", help_heading = "Protein model")]
    idr_local: Option<String>,

    /// Disordered residue ranges excluded from native Gō contacts ([ pairs ]).
    #[arg(long, value_name = "RANGES", help_heading = "Protein model")]
    idr_nonlocal: Option<String>,

    /// RESPAC / user charge override table file (<1-based index> <charge> per line).
    #[arg(long, value_name = "FILE", help_heading = "Protein model")]
    charges: Option<PathBuf>,

    /// Native contact energy scale factor for the Clementi Gō model.
    #[arg(
        long,
        value_name = "FACTOR",
        default_value_t = 1.0,
        help_heading = "Protein model"
    )]
    contact_scale: f64,

    // ---- Geometry ----
    /// Coordinate centering: keep (input coordinates, default), origin, or box.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = ProteinCenterArg::Keep,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    center: ProteinCenterArg,

    /// Vacuum padding in nm added on each side of the periodic box.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    padding: f64,

    // ---- Output files ----
    /// Output directory.
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output files"
    )]
    out: PathBuf,

    /// Base / moleculetype name for output files [default: <stem>_cg].
    #[arg(long, value_name = "NAME", help_heading = "Output files")]
    name: Option<String>,

    /// Write GENESIS .inp control files for min, nvt, and pro.
    #[arg(long, help_heading = "Output files")]
    inp: bool,

    /// Write VMD .vmd visualization scripts (implies --psf).
    #[arg(long, help_heading = "Output files")]
    vmd: bool,

    /// Write a .pdb coordinate file.
    #[arg(long, help_heading = "Output files")]
    pdb: bool,

    /// Write a CHARMM .psf structure file.
    #[arg(long, help_heading = "Output files")]
    psf: bool,

    /// Write a CHARMM .crd coordinate file.
    #[arg(long, help_heading = "Output files")]
    crd: bool,

    /// Write an mmCIF .cif coordinate file.
    #[arg(long, help_heading = "Output files")]
    cif: bool,

    // ---- GENESIS control files (--inp) ----
    /// Simulation temperature in K for GENESIS .inp files [default: 323.15].
    #[arg(
        short = 't',
        long,
        value_name = "K",
        help_heading = "GENESIS control files (--inp)"
    )]
    temperature: Option<f64>,

    /// Minimization (min.inp) overrides as KEY=VAL or <steps> (implies --inp).
    #[arg(
        long,
        value_name = "SPEC",
        help_heading = "GENESIS control files (--inp)"
    )]
    min: Vec<String>,

    /// Equilibration (nvt.inp) overrides as KEY=VAL or <steps> (implies --inp).
    #[arg(
        long,
        value_name = "SPEC",
        help_heading = "GENESIS control files (--inp)"
    )]
    equil: Vec<String>,

    /// Production (pro.inp) overrides as KEY=VAL or <steps> (implies --inp).
    #[arg(
        long,
        value_name = "SPEC",
        help_heading = "GENESIS control files (--inp)"
    )]
    prod: Vec<String>,

    // ---- VMD visualization (--vmd) ----
    /// Use a smooth QuickSurf surface instead of a CA Tube (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    vmd_surface: bool,

    /// Slice the system with a clipping plane (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    vmd_cutaway: bool,

    /// Append a Tachyon ray-tracing render command to each script (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    vmd_render: bool,

    /// VMD background colour name [default: white] (implies --vmd).
    #[arg(long, value_name = "COLOR", help_heading = "VMD visualization (--vmd)")]
    vmd_background: Option<String>,

    /// VMD representation material [default: AOChalky] (implies --vmd).
    #[arg(long, value_name = "NAME", help_heading = "VMD visualization (--vmd)")]
    vmd_material: Option<String>,

    /// Trajectory smoothing window in frames [default: 5] (implies --vmd).
    #[arg(long, value_name = "N", help_heading = "VMD visualization (--vmd)")]
    vmd_smooth: Option<usize>,

    // ---- General ----
    /// Seed value for the RNG (random if unset).
    #[arg(long, value_name = "SEED", help_heading = "General")]
    seed: Option<u64>,

    /// Print only the final summary and any warnings.
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        help_heading = "General"
    )]
    quiet: bool,

    /// Print full file paths in the Write section.
    #[arg(short = 'v', long, help_heading = "General")]
    verbose: bool,

    /// Disable coloured terminal output.
    #[arg(long, help_heading = "General")]
    no_color: bool,

    /// Use ASCII status markers instead of Unicode glyphs.
    #[arg(long, help_heading = "General")]
    ascii: bool,
}

impl ProteinArgs {
    fn want_inp(&self) -> bool {
        self.inp || !self.min.is_empty() || !self.equil.is_empty() || !self.prod.is_empty()
    }

    fn want_vmd(&self) -> bool {
        self.vmd
            || self.vmd_surface
            || self.vmd_cutaway
            || self.vmd_render
            || self.vmd_background.is_some()
            || self.vmd_material.is_some()
            || self.vmd_smooth.is_some()
    }
}

/// Arguments for `isolf lipids`.
#[derive(clap::Args, Debug)]
struct LipidsArgs {
    /// Disable coloured terminal output.
    #[arg(long)]
    no_color: bool,

    /// Use ASCII punctuation instead of Unicode glyphs.
    #[arg(long)]
    ascii: bool,
}

/// Parsed minimization stage settings (before attaching the final box size).
#[derive(Debug, Clone, Copy, PartialEq)]
struct MinSpec {
    num_steps: usize,
    eneout_period: usize,
    crdout_period: usize,
    rstout_period: usize,
    nbupdate_period: usize,
}

impl Default for MinSpec {
    fn default() -> Self {
        Self {
            num_steps: 1000,
            eneout_period: 10,
            crdout_period: 1000,
            rstout_period: 1000,
            nbupdate_period: 10,
        }
    }
}

impl MinSpec {
    fn to_params(self, box_size: [f64; 3], temperature: f64) -> isolf::inp::Minimization {
        isolf::inp::Minimization {
            box_size,
            temperature,
            num_steps: self.num_steps,
            eneout_period: self.eneout_period,
            crdout_period: self.crdout_period,
            rstout_period: self.rstout_period,
            nbupdate_period: self.nbupdate_period,
        }
    }
}

/// Parsed dynamics stage settings (`--equil` or `--prod`).
#[derive(Debug, Clone, Copy, PartialEq)]
struct DynSpec {
    num_steps: usize,
    time_step: f64,
    temperature: Option<f64>,
    eneout_period: usize,
    crdout_period: usize,
    rstout_period: usize,
    stoptr_period: usize,
    nbupdate_period: usize,
    seed: Option<u64>,
}

impl DynSpec {
    fn default_equil() -> Self {
        Self {
            num_steps: 20_000,
            time_step: 0.01,
            temperature: None,
            eneout_period: 100,
            crdout_period: 100,
            rstout_period: 1000,
            stoptr_period: 10,
            nbupdate_period: 10,
            seed: None,
        }
    }

    fn default_prod() -> Self {
        Self {
            num_steps: 500_000,
            time_step: 0.02,
            temperature: None,
            eneout_period: 500,
            crdout_period: 500,
            rstout_period: 1000,
            stoptr_period: 10,
            nbupdate_period: 10,
            seed: None,
        }
    }

    fn to_params(self, default_temperature: f64, seed: u64) -> isolf::inp::Dynamics {
        isolf::inp::Dynamics {
            temperature: self.temperature.unwrap_or(default_temperature),
            num_steps: self.num_steps,
            time_step: self.time_step,
            eneout_period: self.eneout_period,
            crdout_period: self.crdout_period,
            rstout_period: self.rstout_period,
            stoptr_period: self.stoptr_period,
            nbupdate_period: self.nbupdate_period,
            seed: self.seed.unwrap_or(seed),
        }
    }
}

/// Upper bound for a randomly drawn Langevin thermostat seed.
const MAX_ISEED: u64 = 999_999;

/// One optional structure-file format: whether it was requested, its extension,
/// and how to render it from the built system.
type StructureFormat = (bool, &'static str, fn(&EmbeddedSystem) -> String);

/// Colour scheme for the `--help` output.
fn help_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Green.on_default().bold())
        .usage(AnsiColor::Green.on_default().bold())
        .literal(AnsiColor::Cyan.on_default().bold())
        .placeholder(AnsiColor::Cyan.on_default())
}

fn main() {
    // Detected before parsing so `--help` and clap errors also respect `--no-color`.
    let no_color =
        std::env::args().any(|arg| arg == "--no-color") || std::env::var_os("NO_COLOR").is_some();
    let color = if no_color {
        ColorChoice::Never
    } else {
        ColorChoice::Auto
    };
    let command = Cli::command().color(color).styles(help_styles());
    let cli = Cli::from_arg_matches(&command.get_matches()).unwrap_or_else(|e| e.exit());

    let (ascii, quiet, verbose, flag_no_color) = match &cli.command {
        Command::Membrane(a) => (
            a.common.ascii,
            a.common.quiet,
            a.common.verbose,
            a.common.no_color,
        ),
        Command::Vesicle(a) => (
            a.common.ascii,
            a.common.quiet,
            a.common.verbose,
            a.common.no_color,
        ),
        Command::Tubule(a) => (
            a.common.ascii,
            a.common.quiet,
            a.common.verbose,
            a.common.no_color,
        ),
        Command::Protein(a) => (a.ascii, a.quiet, a.verbose, a.no_color),
        Command::Lipids(a) => (a.ascii, false, false, a.no_color),
        Command::Update | Command::Uninstall => (false, false, false, false),
    };
    let color_on = !(no_color || flag_no_color) && std::io::stdout().is_terminal();
    report::init(color_on, !ascii, quiet, verbose);

    if let Err(error) = dispatch(cli.command, color_on, !ascii) {
        report::error(&error.to_string());
        std::process::exit(1);
    }
}

fn dispatch(command: Command, color: bool, unicode: bool) -> Result<(), Box<dyn Error>> {
    match command {
        Command::Membrane(args) => run_membrane(args),
        Command::Vesicle(args) => run_vesicle(args),
        Command::Tubule(args) => run_tubule(args),
        Command::Protein(args) => run_protein(args),
        Command::Lipids(_) => {
            print_lipids(&ForceField::isolf(), color, unicode);
            Ok(())
        }
        Command::Update => update::run(),
        Command::Uninstall => update::uninstall(),
    }
}

/// Validated GENESIS `.inp` stage specs ready to write after the build finishes.
struct ControlSpecs {
    temperature: f64,
    min: MinSpec,
    equil: DynSpec,
    prod: DynSpec,
}

/// Parse and pre-validate `--min`, `--equil`, `--prod`, and `--temperature`
/// before building so bad period multiples fail immediately.
fn prepare_controls(
    common: &CommonBuildArgs,
    geometry: Geometry,
) -> Result<(ControlSpecs, Option<String>), Box<dyn Error>> {
    if let Some(t) = common.temperature
        && (!t.is_finite() || t <= 0.0)
    {
        return Err(format!("--temperature: '{t}' is not a positive temperature in K").into());
    }
    let temperature = common.temperature.unwrap_or(DEFAULT_TEMPERATURE);
    let min = parse_min_spec(&common.min)?;
    let equil = parse_dyn_spec("--equil", &common.equil, DynSpec::default_equil())?;
    let prod = parse_dyn_spec("--prod", &common.prod, DynSpec::default_prod())?;

    let want_inp = common.want_inp();
    if want_inp {
        // Validate step/period divisibility up front (placeholder seed & box).
        isolf::inp::minimization("x", &min.to_params([1.0; 3], temperature))
            .map_err(|e| format!("--min: {e}"))?;
        let check_dyn = |flag: &str, spec: DynSpec| -> Result<(), Box<dyn Error>> {
            let p = spec.to_params(temperature, 1);
            let res = if geometry.is_npt() {
                isolf::inp::npt("x", "x", "x", &p)
            } else {
                isolf::inp::nvt("x", "x", "x", &p)
            };
            res.map(|_| ()).map_err(|e| format!("{flag}: {e}").into())
        };
        check_dyn("--equil", equil)?;
        check_dyn("--prod", prod)?;
    }

    let temp_note = (common.temperature.is_some() && !want_inp).then(|| {
        "--temperature only affects GENESIS .inp files (pass --inp to write them)".to_string()
    });
    Ok((
        ControlSpecs {
            temperature,
            min,
            equil,
            prod,
        },
        temp_note,
    ))
}

fn validate_padding(padding: f64) -> Result<(), Box<dyn Error>> {
    if !padding.is_finite() || padding < 0.0 {
        return Err(format!("--padding: '{padding}' must be a non-negative number of nm").into());
    }
    Ok(())
}

/// Resolve the two leaflet/shell compositions from symmetric `--lipids` or
/// paired (`--upper`/`--lower`, `--outer`/`--inner`) specs.
fn resolve_compositions(
    lipids: &[String],
    first: &[String],
    second: &[String],
    first_flag: &str,
    second_flag: &str,
) -> Result<(Composition, Composition), Box<dyn Error>> {
    let sym = join_specs(lipids);
    if !sym.is_empty() {
        let comp = Composition::parse(&sym)?;
        return Ok((comp.clone(), comp));
    }
    let a = join_specs(first);
    let b = join_specs(second);
    let a_opt = (!a.is_empty()).then_some(a.as_str());
    let b_opt = (!b.is_empty()).then_some(b.as_str());
    let first_spec = a_opt.or(b_opt).ok_or_else(|| {
        format!("missing lipid composition: pass --lipids <SPEC> (or {first_flag}/{second_flag})")
    })?;
    let second_spec = b_opt.or(a_opt).unwrap_or(first_spec);
    Ok((
        Composition::parse(first_spec)?,
        Composition::parse(second_spec)?,
    ))
}

fn phase_notes(common: &CommonBuildArgs, first: &Composition, second: &Composition) -> Vec<String> {
    let mut notes = Vec::new();
    if common.anti_register && common.phase == PhaseArg::Random {
        notes.push(
            "--anti-register has no effect with --phase random (pass --phase stripe or --phase raft)"
                .to_string(),
        );
    } else if common.phase != PhaseArg::Random
        && first.domain_count() <= 1
        && second.domain_count() <= 1
    {
        notes.push(
            "--phase has no effect with a single lipid domain (specify multiple lipids or ';'-separated domains)"
                .to_string(),
        );
    }
    notes
}

/// Parse a finite floating-point number (signed allowed), trimming whitespace.
fn parse_f64(value: &str, what: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("'{value}' is not a valid {what}").into())
}

/// Parse a non-negative finite floating-point number, trimming whitespace.
fn parse_nonneg_f64(value: &str, what: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v >= 0.0)
        .ok_or_else(|| format!("'{value}' is not a non-negative {what}").into())
}

/// Parse repeatable `-p, --protein <PATH>[,KEY=VAL,...]` flags into [`ProteinSpec`]s.
fn parse_protein_specs(specs: &[String]) -> Result<Vec<ProteinSpec>, Box<dyn Error>> {
    let mut out = Vec::with_capacity(specs.len());
    for raw in specs {
        let mut tokens = raw.split(',').map(str::trim).filter(|t| !t.is_empty());
        let path_str = tokens
            .next()
            .ok_or("--protein: missing input structure path")?;
        if path_str.contains('=') {
            return Err(format!(
                "--protein: expected a PDB/mmCIF file path before options, got '{path_str}'"
            )
            .into());
        }

        let mut placement = ProteinPlacement::default();
        let mut options = ProteinOptions::default();
        let mut name_override: Option<String> = None;

        for tok in tokens {
            let (k, v) = tok
                .split_once('=')
                .ok_or_else(|| format!("--protein: expected KEY=VAL option, got '{tok}'"))?;
            let key = k.trim().to_ascii_lowercase();
            let val = v.trim();
            match key.as_str() {
                "z" => {
                    if val.eq_ignore_ascii_case("center") || val.eq_ignore_ascii_case("mid") {
                        placement.z_mode = ZPlacement::Center;
                    } else if val.eq_ignore_ascii_case("keep") {
                        placement.z_mode = ZPlacement::Keep;
                    } else {
                        placement.z_offset = parse_f64(val, "z offset in nm")?;
                    }
                }
                "zoff" | "z_offset" | "dz" => {
                    placement.z_offset = parse_f64(val, "z offset in nm")?;
                }
                "x" => placement.x = Some(parse_f64(val, "x coordinate in nm")?),
                "y" => placement.y = Some(parse_f64(val, "y coordinate in nm")?),
                "tilt" => placement.tilt_deg = parse_f64(val, "tilt angle in degrees")?,
                "rot" | "rotate" => {
                    placement.rot_deg = parse_f64(val, "rotation angle in degrees")?
                }
                "copies" | "n" => placement.copies = parse_pos_usize(val, "copy count")?,
                "clearance" | "margin" => {
                    placement.clearance = parse_nonneg_f64(val, "clearance in nm")?;
                }
                "model" | "ff" => {
                    options.model = match val.to_ascii_lowercase().as_str() {
                        "aicg2+" | "aicg2p" | "aicg2" => ProteinModel::Aicg2Plus,
                        "clementi" | "ccgo" => ProteinModel::ClementiGo,
                        other => {
                            return Err(format!(
                                "--protein: unknown model '{other}' (use aicg2+ or clementi)"
                            )
                            .into());
                        }
                    };
                }
                "scale" => {
                    options.scale = match val.to_ascii_lowercase().as_str() {
                        "general" | "1" => AicgScale::General,
                        "average" | "0" => AicgScale::Average,
                        other => {
                            return Err(format!(
                                "--protein: unknown scale '{other}' (use general or average)"
                            )
                            .into());
                        }
                    };
                }
                "safe-dihedral" | "safe_dihedral" | "safedih" => {
                    options.safe_dihedral = match val.to_ascii_lowercase().as_str() {
                        "0" | "none" => SafeDihedral::None,
                        "1" | "cos2" => SafeDihedral::Cos2,
                        "3" | "sin3" => SafeDihedral::Sin3,
                        other => {
                            return Err(format!(
                                "--protein: unknown safe-dihedral '{other}' (use 0, 1, or 3)"
                            )
                            .into());
                        }
                    };
                }
                "idr" => {
                    let norm = val.replace(['+', ';'], ",");
                    let sel = parse_selection(&norm)?;
                    options.idr_local = sel.clone();
                    options.idr_nonlocal = sel;
                }
                "idr-local" | "idr_local" => {
                    let norm = val.replace(['+', ';'], ",");
                    options.idr_local = parse_selection(&norm)?;
                }
                "idr-nonlocal" | "idr_nonlocal" => {
                    let norm = val.replace(['+', ';'], ",");
                    options.idr_nonlocal = parse_selection(&norm)?;
                }
                "charges" | "respac" => {
                    let content = fs::read_to_string(val).map_err(|e| {
                        format!("--protein: failed to read charges file '{val}': {e}")
                    })?;
                    options.charges = parse_charge_table(&content)?;
                }
                "contact-scale" | "contact_scale" => {
                    options.ccgo_contact_scale = parse_pos_f64(val, "contact scale")?;
                }
                "name" => name_override = Some(val.to_string()),
                other => {
                    return Err(format!(
                        "--protein: unknown key '{other}' (use z, x, y, tilt, rot, copies, clearance, model, scale, safe-dihedral, idr, idr-local, idr-nonlocal, charges, name)"
                    )
                    .into());
                }
            }
        }

        let protein = Protein::from_file(Path::new(path_str), name_override.as_deref(), &options)?;
        out.push(ProteinSpec { protein, placement });
    }
    Ok(out)
}

fn run_membrane(args: MembraneArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    validate_padding(args.padding)?;

    // Parse & validate all arguments before printing the banner.
    let (upper, lower) =
        resolve_compositions(&args.lipids, &args.upper, &args.lower, "--upper", "--lower")?;
    let counts = parse_count_spec(&args.count)?;
    let box_xy = parse_box_spec(&args.box_size)?;
    if counts.is_none() && box_xy.is_none() {
        return Err("missing membrane size: pass --count <N>, --box <L>, or both".into());
    }
    let (sizing, balance_note) = resolve_membrane_sizing(
        args.no_balance,
        &force_field,
        &upper,
        &lower,
        counts,
        box_xy,
    )
    .map_err(|e| enrich_build_error(&force_field, e))?;
    let protein_specs = parse_protein_specs(&args.common.proteins)?;
    let (controls, temp_note) = prepare_controls(&args.common, Geometry::Membrane)?;

    let seed = args.common.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);
    let options = BuildOptions {
        padding: args.padding,
        phase: args.common.phase.into(),
        anti_register: args.common.anti_register,
    };

    report::header();
    report::section("Build");
    let spin = report::spin("building membrane");
    let built = if args.common.want_gif() {
        embed_membrane_recorded(
            &force_field,
            Geometry::Membrane.title(),
            &upper,
            &lower,
            sizing,
            &options,
            &protein_specs,
            &mut rng,
        )
    } else {
        embed_membrane(
            &force_field,
            Geometry::Membrane.title(),
            &upper,
            &lower,
            sizing,
            &options,
            &protein_specs,
            &mut rng,
        )
        .map(|s| (s, Vec::new()))
    };
    spin.finish();
    let (system, gif_frames) = built.map_err(|e| enrich_build_error(&force_field, e))?;

    let mut notes = phase_notes(&args.common, &upper, &lower);
    if let Some(n) = balance_note {
        notes.push(n);
    }
    if let Some(n) = temp_note {
        notes.push(n);
    }
    for spec in &protein_specs {
        notes.extend(spec.protein.warnings.iter().cloned());
    }
    notes.extend(system.warnings.iter().cloned());
    let gif_bounds = [0.0, 0.0, system.box_size[0], system.box_size[1]];
    finish_build(
        system,
        gif_frames,
        gif_bounds,
        &force_field,
        &args.common,
        Geometry::Membrane,
        args.center,
        &controls,
        &notes,
        seed,
        &mut rng,
    )
}

fn run_vesicle(args: VesicleArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    validate_padding(args.padding)?;

    // Parse & validate all arguments before printing the banner.
    let (outer, inner) =
        resolve_compositions(&args.lipids, &args.outer, &args.inner, "--outer", "--inner")?;
    let radius = parse_radius_spec(&args.radius)?
        .ok_or("missing vesicle radius: pass --radius <R> (e.g. 20, ri=15, or ri=14,ro=20)")?;
    // Validate lipid names before printing the banner.
    isolf::membrane::area_per_lipid(&force_field, &outer)
        .map_err(|e| enrich_build_error(&force_field, e))?;
    isolf::membrane::area_per_lipid(&force_field, &inner)
        .map_err(|e| enrich_build_error(&force_field, e))?;
    let protein_specs = parse_protein_specs(&args.common.proteins)?;
    let (controls, temp_note) = prepare_controls(&args.common, Geometry::Vesicle)?;

    let seed = args.common.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);
    let options = BuildOptions {
        padding: args.padding,
        phase: args.common.phase.into(),
        anti_register: args.common.anti_register,
    };

    report::header();
    report::section("Build");
    let spin = report::spin("building vesicle");
    let built = if args.common.want_gif() {
        embed_vesicle_recorded(
            &force_field,
            Geometry::Vesicle.title(),
            &outer,
            &inner,
            radius,
            &options,
            &protein_specs,
            &mut rng,
        )
    } else {
        embed_vesicle(
            &force_field,
            Geometry::Vesicle.title(),
            &outer,
            &inner,
            radius,
            &options,
            &protein_specs,
            &mut rng,
        )
        .map(|s| (s, Vec::new()))
    };
    spin.finish();
    let (system, gif_frames) = built.map_err(|e| enrich_build_error(&force_field, e))?;

    let mut notes = phase_notes(&args.common, &outer, &inner);
    if let Some(n) = temp_note {
        notes.push(n);
    }
    for spec in &protein_specs {
        notes.extend(spec.protein.warnings.iter().cloned());
    }
    notes.extend(system.warnings.iter().cloned());
    finish_build(
        system,
        gif_frames,
        [-2.0, -2.0, 2.0, 2.0],
        &force_field,
        &args.common,
        Geometry::Vesicle,
        args.center,
        &controls,
        &notes,
        seed,
        &mut rng,
    )
}

fn run_tubule(args: TubuleArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    validate_padding(args.padding)?;
    if !args.length.is_finite() || args.length <= 0.0 {
        return Err(format!(
            "--length: '{}' must be a positive number of nm",
            args.length
        )
        .into());
    }

    let (outer, inner) =
        resolve_compositions(&args.lipids, &args.outer, &args.inner, "--outer", "--inner")?;
    let radius = match parse_radius_spec(&args.radius)?
        .ok_or("missing tubule radius: pass --radius <R> (e.g. 10, ri=6, or ri=6,ro=10.5)")?
    {
        VesicleRadius::Inner(r) => TubuleRadius::Inner(r),
        VesicleRadius::Outer(r) => TubuleRadius::Outer(r),
        VesicleRadius::Both { inner, outer } => TubuleRadius::Both { inner, outer },
    };
    isolf::membrane::area_per_lipid(&force_field, &outer)
        .map_err(|e| enrich_build_error(&force_field, e))?;
    isolf::membrane::area_per_lipid(&force_field, &inner)
        .map_err(|e| enrich_build_error(&force_field, e))?;
    let protein_specs = parse_protein_specs(&args.common.proteins)?;
    let (controls, temp_note) = prepare_controls(&args.common, Geometry::Tubule)?;

    let seed = args.common.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);
    let options = BuildOptions {
        padding: args.padding,
        phase: args.common.phase.into(),
        anti_register: args.common.anti_register,
    };

    report::header();
    report::section("Build");
    let spin = report::spin("building tubule");
    let built = if args.common.want_gif() {
        embed_tubule_recorded(
            &force_field,
            Geometry::Tubule.title(),
            &outer,
            &inner,
            radius,
            args.length,
            &options,
            &protein_specs,
            &mut rng,
        )
    } else {
        embed_tubule(
            &force_field,
            Geometry::Tubule.title(),
            &outer,
            &inner,
            radius,
            args.length,
            &options,
            &protein_specs,
            &mut rng,
        )
        .map(|s| (s, Vec::new()))
    };
    spin.finish();
    let (system, gif_frames) = built.map_err(|e| enrich_build_error(&force_field, e))?;

    let mut notes = phase_notes(&args.common, &outer, &inner);
    if let Some(n) = temp_note {
        notes.push(n);
    }
    for spec in &protein_specs {
        notes.extend(spec.protein.warnings.iter().cloned());
    }
    notes.extend(system.warnings.iter().cloned());
    let outer_radius = (0.5 * system.box_size[0] - args.padding).max(1e-6);
    let gif_bounds = [0.0, 0.0, TAU * outer_radius, system.box_size[2]];
    finish_build(
        system,
        gif_frames,
        gif_bounds,
        &force_field,
        &args.common,
        Geometry::Tubule,
        args.center,
        &controls,
        &notes,
        seed,
        &mut rng,
    )
}

/// Standalone `isolf protein <INPUT>` coarse-graining workflow (AICG2+ / Clementi Gō).
fn run_protein(args: ProteinArgs) -> Result<(), Box<dyn Error>> {
    validate_padding(args.padding)?;
    if !args.contact_scale.is_finite() || args.contact_scale <= 0.0 {
        return Err(format!(
            "--contact-scale: '{}' must be a positive number",
            args.contact_scale
        )
        .into());
    }
    if let Some(t) = args.temperature
        && (!t.is_finite() || t <= 0.0)
    {
        return Err(format!("--temperature: '{t}' is not a positive temperature in K").into());
    }

    let mut idr_local = match &args.idr {
        Some(s) => parse_selection(s)?,
        None => Vec::new(),
    };
    let mut idr_nonlocal = idr_local.clone();
    if let Some(s) = &args.idr_local {
        idr_local.extend(parse_selection(s)?);
        idr_local.sort_unstable();
        idr_local.dedup();
    }
    if let Some(s) = &args.idr_nonlocal {
        idr_nonlocal.extend(parse_selection(s)?);
        idr_nonlocal.sort_unstable();
        idr_nonlocal.dedup();
    }
    let charges = match &args.charges {
        Some(path) => {
            let content = fs::read_to_string(path)
                .map_err(|e| format!("--charges: failed to read {}: {e}", path.display()))?;
            parse_charge_table(&content)?
        }
        None => Vec::new(),
    };

    let temperature = args.temperature.unwrap_or(DEFAULT_TEMPERATURE);
    let min = parse_min_spec(&args.min)?;
    let equil = parse_dyn_spec("--equil", &args.equil, DynSpec::default_equil())?;
    let prod = parse_dyn_spec("--prod", &args.prod, DynSpec::default_prod())?;
    let want_inp = args.want_inp();
    let want_vmd = args.want_vmd();
    if want_inp {
        isolf::inp::minimization("x", &min.to_params([1.0; 3], temperature))
            .map_err(|e| format!("--min: {e}"))?;
        isolf::inp::nvt("x", "x", "x", &equil.to_params(temperature, 1))
            .map_err(|e| format!("--equil: {e}"))?;
        isolf::inp::nvt("x", "x", "x", &prod.to_params(temperature, 1))
            .map_err(|e| format!("--prod: {e}"))?;
    }

    let prot_options = ProteinOptions {
        model: args.model.into(),
        scale: args.scale.into(),
        safe_dihedral: args.safe_dihedral.into(),
        idr_local,
        idr_nonlocal,
        charges,
        ccgo_contact_scale: args.contact_scale,
    };

    let seed = args.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);

    report::header();
    report::section("Build");
    let spin = report::spin("coarse-graining protein");
    let mut protein = Protein::from_file(&args.input, args.name.as_deref(), &prot_options)?;
    spin.finish();

    let box_size = protein.box_with_padding(args.padding);
    match args.center {
        ProteinCenterArg::Keep => {}
        ProteinCenterArg::Origin => {
            let c = protein.centroid();
            protein.translate([-c[0], -c[1], -c[2]]);
        }
        ProteinCenterArg::Box => {
            let c = protein.centroid();
            protein.translate([
                0.5 * box_size[0] - c[0],
                0.5 * box_size[1] - c[1],
                0.5 * box_size[2] - c[2],
            ]);
        }
    }

    let model_tag = match protein.model {
        ProteinModel::Aicg2Plus => "AICG2+",
        ProteinModel::ClementiGo => "Clementi Gō",
    };
    let chain_word = if protein.chains.len() == 1 {
        "chain"
    } else {
        "chains"
    };
    report::ok(
        "protein",
        &format!(
            "{} ({}) {} {} {chain_word} {} {} beads {} {} contacts {} Q = {:+.1} e",
            protein.name,
            model_tag,
            report::dot(),
            protein.chains.len(),
            report::dot(),
            report::thousands(protein.beads.len()),
            report::dot(),
            report::thousands(protein.pairs.len()),
            report::dot(),
            protein.net_charge(),
        ),
    );
    for w in &protein.warnings {
        report::warn(w);
    }
    if args.temperature.is_some() && !want_inp {
        report::warn("--temperature only affects GENESIS .inp files (pass --inp to write them)");
    }

    fs::create_dir_all(&args.out)?;
    let base = protein.name.clone();
    let named = |ext: &str| args.out.join(format!("{base}.{ext}"));

    write(&named("gro"), &protein.to_gro(box_size))?;
    write(&named("itp"), &protein.to_itp())?;
    write(&args.out.join("aicg2p.itp"), &render_aicg2p_itp())?;
    write(&named("top"), &protein.to_top())?;
    if args.pdb {
        write(&named("pdb"), &protein.to_pdb())?;
    }
    if args.psf || want_vmd {
        write(&named("psf"), &protein.to_psf())?;
    }
    if args.crd {
        write(&named("crd"), &protein.to_crd())?;
    }
    if args.cif {
        write(&named("cif"), &protein.to_cif())?;
    }

    let controls = ControlSpecs {
        temperature,
        min,
        equil,
        prod,
    };
    if want_inp {
        write_control_files(
            &args.out,
            &base,
            Geometry::Vesicle,
            box_size,
            &controls,
            &mut rng,
        )?;
    }
    if want_vmd {
        let ff = ForceField::isolf();
        let particles = protein.to_particles(1, 1);
        let background = args.vmd_background.as_deref().unwrap_or("white");
        let material = args.vmd_material.as_deref().unwrap_or("AOChalky");
        let vmd_opts = VmdOptions {
            surface: args.vmd_surface,
            cutaway: args.vmd_cutaway,
            render: args.vmd_render,
            background,
            material,
            smooth: args.vmd_smooth.unwrap_or(5),
            lipid_coloring: isolf::vmd::LipidColoring::Head,
            recenter_to_box: matches!(args.center, ProteinCenterArg::Origin),
        };
        let script = |source| {
            isolf::vmd::vmd_script(&ff, &particles, &base, Layout::Protein, source, &vmd_opts)
        };
        write(
            &args.out.join(format!("{base}.vmd")),
            &script(VmdSource::Structure),
        )?;
        if want_inp {
            for (stage, steps, crdout) in [
                ("min", controls.min.num_steps, controls.min.crdout_period),
                (
                    "nvt",
                    controls.equil.num_steps,
                    controls.equil.crdout_period,
                ),
                ("pro", controls.prod.num_steps, controls.prod.crdout_period),
            ] {
                let frames = (crdout > 0).then(|| steps / crdout);
                write(
                    &args.out.join(format!("{stage}.vmd")),
                    &script(VmdSource::Trajectory {
                        name: stage,
                        frames,
                    }),
                )?;
            }
        }
    }

    let log_path = named("log");
    let log_content = render_protein_log(&protein, box_size, seed, args.seed.is_some());
    write(&log_path, &log_content)?;

    let out_dir = args.out.display().to_string();
    report::write_section(&out_dir);
    let next = want_inp.then(|| {
        let cd = if out_dir == "." {
            String::new()
        } else {
            format!("cd {out_dir} && ")
        };
        format!("run in GENESIS: {cd}atdyn min.inp (then nvt.inp, pro.inp)")
    });
    let prot_summary = [(protein.name.clone(), 1usize)];
    report::summary(&report::Summary {
        identity: "protein",
        box_nm: box_size,
        proteins: &prot_summary,
        lipids: &[],
        leaflets: None,
        beads: protein.beads.len(),
        next,
    });
    Ok(())
}

/// Format a detailed protein coarse-graining log (matching `genesis_cg_tool`'s
/// chain table, Rg/Rc, net charge, contact count, and reproducible invocation).
fn render_protein_log(
    protein: &Protein,
    box_size: [f64; 3],
    seed: u64,
    seed_given: bool,
) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "command  {}", reproduce_command(seed, seed_given));
    let _ = writeln!(
        s,
        "box      {:.4} x {:.4} x {:.4} nm",
        box_size[0], box_size[1], box_size[2]
    );
    let _ = writeln!(s);
    let _ = writeln!(s, "Input structure (atomic):");
    let _ = writeln!(s, "  atoms      : {}", protein.source_atoms);
    let _ = writeln!(s, "  residues   : {}", protein.source_residues);
    let _ = writeln!(s, "  chains     : {}", protein.chains.len());
    let _ = writeln!(s);
    let _ = writeln!(s, "Chain info (CG):");
    let _ = writeln!(
        s,
        "  | Chain | ID | # bead | start --   end |  Rg (nm) | net charge |"
    );
    let _ = writeln!(
        s,
        "  |-------+----+--------+----------------+----------+------------|"
    );
    for (ci, ch) in protein.chains.iter().enumerate() {
        let cid = if ch.chain_id == ' ' { '-' } else { ch.chain_id };
        let _ = writeln!(
            s,
            "  | {:>5} | {:>2} | {:>6} | {:>5} -- {:>5} | {:>8.3} | {:>+10.3} |",
            ci + 1,
            cid,
            ch.range.len(),
            ch.range.start + 1,
            ch.range.end,
            protein.chain_rg(ci),
            protein.chain_charge(ci),
        );
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "CG model summary:");
    let _ = writeln!(s, "  beads      : {}", protein.beads.len());
    let _ = writeln!(s, "  bonds      : {}", protein.bonds.len());
    let _ = writeln!(s, "  angles     : {}", protein.angles.len());
    let _ = writeln!(s, "  dihedrals  : {}", protein.dihedrals.len());
    let _ = writeln!(s, "  contacts   : {}", protein.pairs.len());
    let _ = writeln!(s, "  Rg         : {:.3} nm", protein.radius_of_gyration());
    let _ = writeln!(s, "  Rc         : {:.3} nm", protein.circumsphere_radius());
    let _ = writeln!(s, "  net charge : {:+.3} e", protein.net_charge());
    s
}

#[allow(clippy::too_many_arguments)]
fn finish_build(
    mut system: EmbeddedSystem,
    gif_frames: Vec<RelaxFrame>,
    gif_bounds: [f64; 4],
    force_field: &ForceField,
    common: &CommonBuildArgs,
    geometry: Geometry,
    center: CenterArg,
    controls: &ControlSpecs,
    notes: &[String],
    seed: u64,
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    let identity = geometry.label();
    let total_removed = system.carved_lipids + system.pruned_lipids;
    let detail = if system.has_proteins() {
        let prot_word = if system.total_proteins() == 1 {
            "protein"
        } else {
            "proteins"
        };
        let carved_str = if total_removed > 0 {
            format!(" ({} carved)", report::thousands(total_removed))
        } else {
            String::new()
        };
        format!(
            "{} lipids{carved_str} {} {} {prot_word} {} {} beads",
            report::thousands(system.total_lipids()),
            report::dot(),
            report::thousands(system.total_proteins()),
            report::dot(),
            report::thousands(system.particles.len()),
        )
    } else {
        format!(
            "{} lipids {} {} beads",
            report::thousands(system.total_lipids()),
            report::dot(),
            report::thousands(system.particles.len()),
        )
    };
    report::ok(identity, &detail);
    for note in notes {
        report::warn(note);
    }

    let base = common
        .name
        .clone()
        .unwrap_or_else(|| geometry.label().to_string());

    write_outputs(
        &mut system,
        force_field,
        common,
        geometry,
        center,
        controls,
        &base,
        &gif_frames,
        gif_bounds,
        rng,
    )?;

    let out_dir = common.out.display().to_string();
    let box_size = system.box_size;
    let beads = system.particles.len();

    let files: Vec<String> = report::recorded_files()
        .iter()
        .filter_map(|p| p.file_name().map(|f| f.to_string_lossy().into_owned()))
        .collect();
    let command = reproduce_command(seed, common.seed.is_some());
    let log_path = common.out.join(format!("{base}.log"));
    runlog::Run {
        command: &command,
        seed,
        duration_secs: report::elapsed_secs(),
        identity,
        box_nm: box_size,
        beads,
        out_dir: &out_dir,
        files: &files,
    }
    .write(&log_path)?;
    report::record(&log_path);

    report::write_section(&out_dir);
    let next = common.want_inp().then(|| {
        let equil = if geometry.is_npt() {
            "npt.inp"
        } else {
            "nvt.inp"
        };
        let cd = if out_dir == "." {
            String::new()
        } else {
            format!("cd {out_dir} && ")
        };
        format!("run in GENESIS: {cd}atdyn min.inp (then {equil}, pro.inp)")
    });
    let leaflets = leaflet_split(&system, geometry);
    report::summary(&report::Summary {
        identity,
        box_nm: box_size,
        proteins: system.protein_counts.as_slice(),
        lipids: system.lipid_counts.as_slice(),
        leaflets: Some(&leaflets),
        beads,
        next,
    });
    Ok(())
}

/// Turn an [`UnknownLipid`](BuildError::UnknownLipid) into an actionable message
/// with close matches and a pointer to `isolf lipids`.
fn enrich_build_error(force_field: &ForceField, err: BuildError) -> Box<dyn Error> {
    if let BuildError::UnknownLipid(ref name) = err {
        let suggestions = suggest_lipids(force_field, name);
        if suggestions.is_empty() {
            return format!(
                "unknown lipid '{name}' (run 'isolf lipids' to list all {} supported lipids)",
                force_field.lipids.len()
            )
            .into();
        }
        return format!(
            "unknown lipid '{name}' — did you mean {}? (run 'isolf lipids' to list all {})",
            suggestions.join(", "),
            force_field.lipids.len()
        )
        .into();
    }
    err.into()
}

/// Suggest up to 5 closest lipid names within Levenshtein edit distance 2.
fn suggest_lipids<'a>(force_field: &'a ForceField, query: &str) -> Vec<&'a str> {
    let q = query.to_ascii_uppercase();
    let mut ranked: Vec<(usize, &str)> = force_field
        .lipids
        .iter()
        .filter_map(|l| {
            let d = edit_distance(&q, &l.name);
            (d <= 2).then_some((d, l.name.as_str()))
        })
        .collect();
    ranked.sort_by_key(|&(d, name)| (d, name));
    if let Some(&(best, _)) = ranked.first() {
        ranked.retain(|&(d, _)| d == best);
    }
    ranked.into_iter().take(5).map(|(_, n)| n).collect()
}

fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut curr = vec![0; b.len() + 1];
    for (i, &ca) in a.iter().enumerate() {
        curr[0] = i + 1;
        for (j, &cb) in b.iter().enumerate() {
            let cost = usize::from(ca != cb);
            curr[j + 1] = (prev[j + 1] + 1).min(curr[j] + 1).min(prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[b.len()]
}

/// Print the 7×5 iSoLF lipid table directly from the force field.
fn print_lipids(ff: &ForceField, color: bool, unicode: bool) {
    let paint = |s: &str, code: &str| {
        if color {
            format!("\x1b[{code}m{s}\x1b[0m")
        } else {
            s.to_string()
        }
    };
    let tails = [
        ("DL", "lauroyl (12:0)"),
        ("DM", "myristoyl (14:0)"),
        ("DP", "palmitoyl (16:0)"),
        ("DS", "stearoyl (18:0)"),
        ("DO", "dioleoyl (18:1)"),
        ("PO", "16:0 / 18:1"),
        ("SO", "18:0 / 18:1"),
    ];
    let heads = ["PC", "PE", "PG", "PS", "PA"];
    let sep = if unicode { "×" } else { "x" };

    println!(
        "{} ({} species: {} {sep} {}, named <tail><head>)\n",
        paint("iSoLF lipids", "1;36"),
        ff.lipids.len(),
        tails.len(),
        heads.len(),
    );

    // Header row with each head group's net charge computed from the force field.
    let head_cols: Vec<String> = heads
        .iter()
        .map(|h| {
            let sample = format!("PO{h}");
            let q: f64 = ff
                .lipid(&sample)
                .map(|l| {
                    l.beads
                        .iter()
                        .filter_map(|b| ff.bead_type(b))
                        .map(|b| b.charge)
                        .sum()
                })
                .unwrap_or(0.0);
            let q_label = if q.abs() < 0.01 {
                "0".to_string()
            } else {
                format!("{q:+.0}")
            };
            format!("{h} ({q_label})")
        })
        .collect();

    println!(
        "  {}  {}  {}  {}  {}  {}  {}  {}",
        paint(&format!("{:<4}", "Tail"), "1"),
        paint(&format!("{:<18}", "Acyl chains"), "1"),
        paint(&format!("{:<8}", head_cols[0]), "1"),
        paint(&format!("{:<8}", head_cols[1]), "1"),
        paint(&format!("{:<8}", head_cols[2]), "1"),
        paint(&format!("{:<8}", head_cols[3]), "1"),
        paint(&format!("{:<8}", head_cols[4]), "1"),
        paint("Beads", "1"),
    );

    for (tail, desc) in tails {
        let names: Vec<String> = heads.iter().map(|h| format!("{tail}{h}")).collect();
        let pc_beads = ff.lipid(&names[0]).map_or(0, |l| l.beads.len());
        let pa_beads = ff.lipid(&names[4]).map_or(0, |l| l.beads.len());
        println!(
            "  {}  {:<18}  {:<8}  {:<8}  {:<8}  {:<8}  {:<8}  {}",
            paint(&format!("{tail:<4}"), "36"),
            desc,
            names[0],
            names[1],
            names[2],
            names[3],
            names[4],
            paint(&format!("{pc_beads} (PA: {pa_beads})"), "2"),
        );
    }
    println!(
        "\n  {}  isolf membrane -L \"POPC=3,POPS=1\" -n 256 --top -o my_membrane",
        paint("Example:", "2")
    );
}

/// Reconstruct the invocation as a re-runnable, shell-quoted line, with the seed
/// pinned in (appended when the user did not pass `--seed`).
fn reproduce_command(seed: u64, seed_given: bool) -> String {
    let mut parts: Vec<String> = std::env::args().map(|a| shell_quote(&a)).collect();
    if let Some(first) = parts.first_mut() {
        *first = "isolf".to_string(); // normalise argv[0], which may be a path
    }
    if !seed_given {
        parts.push("--seed".to_string());
        parts.push(seed.to_string());
    }
    parts.join(" ")
}

/// Single-quote an argument for a POSIX shell when it holds whitespace or shell
/// metacharacters, so the logged command pastes back verbatim.
fn shell_quote(arg: &str) -> String {
    let needs_quote = arg.is_empty()
        || arg
            .chars()
            .any(|c| c.is_whitespace() || "\"'\\$`*?[]{}()<>|&;#~!".contains(c));
    if needs_quote {
        format!("'{}'", arg.replace('\'', r"'\''"))
    } else {
        arg.to_string()
    }
}

/// What is being built, which selects the production control file, the camera,
/// and the default output naming.
#[derive(Clone, Copy, PartialEq)]
enum Geometry {
    Membrane,
    Vesicle,
    Tubule,
}

impl Geometry {
    /// Short mode label (`"membrane"`, `"vesicle"`, or `"tubule"`), also used as
    /// the default output base name when `--name` is not given.
    fn label(self) -> &'static str {
        match self {
            Geometry::Membrane => "membrane",
            Geometry::Vesicle => "vesicle",
            Geometry::Tubule => "tubule",
        }
    }

    /// System title written into the structure-file headers.
    fn title(self) -> &'static str {
        match self {
            Geometry::Membrane => DEFAULT_NAME,
            Geometry::Vesicle => "CG vesicle model",
            Geometry::Tubule => "CG tubule model",
        }
    }

    /// VMD camera layout (curved geometries get perspective views).
    fn layout(self) -> Layout {
        match self {
            Geometry::Membrane => Layout::Membrane,
            Geometry::Vesicle => Layout::Vesicle,
            Geometry::Tubule => Layout::Tubule,
        }
    }

    /// Whether the production/equilibration ensemble is NPT (a flat membrane) as
    /// opposed to NVT (a vesicle or tubule, at fixed volume).
    fn is_npt(self) -> bool {
        self == Geometry::Membrane
    }
}

/// Parse a single positive finite floating-point number, trimming whitespace.
fn parse_pos_f64(value: &str, unit: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v > 0.0)
        .ok_or_else(|| format!("'{value}' is not a positive {unit}").into())
}

/// Parse a positive integer (`>= 1`), trimming whitespace.
fn parse_pos_usize(value: &str, what: &str) -> Result<usize, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<usize>()
        .ok()
        .filter(|n| *n > 0)
        .ok_or_else(|| format!("'{value}' is not a positive {what}").into())
}

/// Parse a non-negative integer (`>= 0`), trimming whitespace.
fn parse_usize(value: &str, what: &str) -> Result<usize, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<usize>()
        .map_err(|_| format!("'{value}' is not a non-negative {what}").into())
}

/// Walk non-empty comma-separated `KEY=VALUE` (or bare `VALUE`) tokens across
/// repeatable CLI flag occurrences. Returns `true` if any token was seen.
fn for_each_kv(
    specs: &[String],
    mut handle: impl FnMut(Option<&str>, &str) -> Result<(), Box<dyn Error>>,
) -> Result<bool, Box<dyn Error>> {
    let mut seen = false;
    for token in specs.iter().flat_map(|s| s.split(',')) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        seen = true;
        match token.split_once('=') {
            Some((k, v)) => handle(Some(k.trim()), v)?,
            None => handle(None, token)?,
        }
    }
    Ok(seen)
}

/// Parse `--count` occurrences into `(upper, lower)` lipid counts.
fn parse_count_spec(specs: &[String]) -> Result<Option<(usize, usize)>, Box<dyn Error>> {
    let (mut up, mut lo) = (None, None);
    let seen = for_each_kv(specs, |key, val| match key {
        None => {
            up = Some(parse_pos_usize(val, "lipid count")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("up") || k.eq_ignore_ascii_case("upper") => {
            up = Some(parse_pos_usize(val, "lipid count")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("lo") || k.eq_ignore_ascii_case("lower") => {
            lo = Some(parse_pos_usize(val, "lipid count")?);
            Ok(())
        }
        Some(other) => Err(format!("--count: unknown key '{other}' (use up or lo)").into()),
    })?;
    if !seen {
        return Ok(None);
    }
    match (up, lo) {
        (Some(u), Some(l)) => Ok(Some((u, l))),
        (Some(n), None) | (None, Some(n)) => Ok(Some((n, n))),
        (None, None) => Err("--count needs a value, e.g. `256` or `up=300,lo=200`".into()),
    }
}

/// Parse `--box` occurrences into `(x, y)` dimensions in nm.
fn parse_box_spec(specs: &[String]) -> Result<Option<(f64, f64)>, Box<dyn Error>> {
    let (mut x, mut y) = (None, None);
    let seen = for_each_kv(specs, |key, val| match key {
        None => {
            x = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("x") => {
            x = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("y") => {
            y = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(other) => Err(format!("--box: unknown key '{other}' (use x or y)").into()),
    })?;
    if !seen {
        return Ok(None);
    }
    match (x, y) {
        (Some(x), Some(y)) => Ok(Some((x, y))),
        (Some(s), None) | (None, Some(s)) => Ok(Some((s, s))),
        (None, None) => Err("--box needs a value, e.g. `20` or `x=30,y=15`".into()),
    }
}

/// Parse `--radius` occurrences into a [`VesicleRadius`].
fn parse_radius_spec(specs: &[String]) -> Result<Option<VesicleRadius>, Box<dyn Error>> {
    let (mut inner, mut outer) = (None, None);
    let seen = for_each_kv(specs, |key, val| match key {
        None => {
            outer = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("ri") || k.eq_ignore_ascii_case("inner") => {
            inner = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(k) if k.eq_ignore_ascii_case("ro") || k.eq_ignore_ascii_case("outer") => {
            outer = Some(parse_pos_f64(val, "number of nm")?);
            Ok(())
        }
        Some(other) => Err(format!("--radius: unknown key '{other}' (use ro or ri)").into()),
    })?;
    if !seen {
        return Ok(None);
    }
    match (inner, outer) {
        (Some(inner), Some(outer)) => Ok(Some(VesicleRadius::Both { inner, outer })),
        (Some(inner), None) => Ok(Some(VesicleRadius::Inner(inner))),
        (None, Some(outer)) => Ok(Some(VesicleRadius::Outer(outer))),
        (None, None) => Err("--radius needs a value, e.g. `20` or `ri=14,ro=20`".into()),
    }
}

/// Parse `--min` override tokens into a [`MinSpec`].
fn parse_min_spec(specs: &[String]) -> Result<MinSpec, Box<dyn Error>> {
    let mut out = MinSpec::default();
    for_each_kv(specs, |key, val| match key {
        None => {
            out.num_steps = parse_pos_usize(val, "step count")?;
            Ok(())
        }
        Some(k) => match k.to_ascii_lowercase().as_str() {
            "steps" | "nsteps" | "num_steps" => {
                out.num_steps = parse_pos_usize(val, "step count")?;
                Ok(())
            }
            "eneout" | "eneout_period" => {
                out.eneout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "crdout" | "crdout_period" => {
                out.crdout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "rstout" | "rstout_period" => {
                out.rstout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "nbupdate" | "nbupdate_period" => {
                out.nbupdate_period = parse_usize(val, "period")?;
                Ok(())
            }
            other => Err(format!(
                "--min: unknown key '{other}' (use steps, eneout, crdout, rstout, nbupdate)"
            )
            .into()),
        },
    })?;
    Ok(out)
}

/// Parse `--equil` or `--prod` override tokens into a [`DynSpec`].
fn parse_dyn_spec(
    flag: &str,
    specs: &[String],
    mut out: DynSpec,
) -> Result<DynSpec, Box<dyn Error>> {
    for_each_kv(specs, |key, val| {
        match key {
        None => {
            out.num_steps = parse_pos_usize(val, "step count")?;
            Ok(())
        }
        Some(k) => match k.to_ascii_lowercase().as_str() {
            "steps" | "nsteps" | "num_steps" => {
                out.num_steps = parse_pos_usize(val, "step count")?;
                Ok(())
            }
            "dt" | "timestep" | "time_step" => {
                out.time_step = parse_pos_f64(val, "time step")?;
                Ok(())
            }
            "temp" | "temperature" => {
                out.temperature = Some(parse_pos_f64(val, "temperature in K")?);
                Ok(())
            }
            "eneout" | "eneout_period" => {
                out.eneout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "crdout" | "crdout_period" => {
                out.crdout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "rstout" | "rstout_period" => {
                out.rstout_period = parse_usize(val, "period")?;
                Ok(())
            }
            "stoptr" | "stoptr_period" => {
                out.stoptr_period = parse_usize(val, "period")?;
                Ok(())
            }
            "nbupdate" | "nbupdate_period" => {
                out.nbupdate_period = parse_usize(val, "period")?;
                Ok(())
            }
            "seed" | "iseed" => {
                let s = val
                    .trim()
                    .parse::<u64>()
                    .map_err(|_| format!("'{val}' is not a valid seed"))?;
                out.seed = Some(s);
                Ok(())
            }
            other => Err(format!(
                "{flag}: unknown key '{other}' (use steps, dt, temp, eneout, crdout, rstout, stoptr, nbupdate, seed)"
            )
            .into()),
        },
    }
    })?;
    Ok(out)
}

/// Combine repeatable, comma-separated composition flags into one string.
fn join_specs(specs: &[String]) -> String {
    specs
        .iter()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join(",")
}

/// Resolve the flat-membrane sizing, matching the two leaflets' areas in count
/// mode unless `no_balance` is set. Also validates all lipid names up front.
fn resolve_membrane_sizing(
    no_balance: bool,
    force_field: &ForceField,
    upper: &Composition,
    lower: &Composition,
    counts: Option<(usize, usize)>,
    box_xy: Option<(f64, f64)>,
) -> Result<(Sizing, Option<String>), BuildError> {
    let upper_apl = isolf::membrane::area_per_lipid(force_field, upper)?;
    let lower_apl = isolf::membrane::area_per_lipid(force_field, lower)?;
    match (counts, box_xy) {
        (Some((upper_count, lower_count)), None) => {
            let (u, l, note) =
                balanced_counts(upper_apl, lower_apl, upper_count, lower_count, !no_balance);
            Ok((Sizing::Count { upper: u, lower: l }, note))
        }
        (None, Some((x, y))) => Ok((Sizing::Box { x, y }, None)),
        (Some((upper, lower)), Some((x, y))) => {
            Ok((Sizing::CountInBox { upper, lower, x, y }, None))
        }
        (None, None) => unreachable!("membrane requires --count or --box"),
    }
}

/// Resolve per-leaflet counts so both leaflets start at matched area.
fn balanced_counts(
    upper_apl: f64,
    lower_apl: f64,
    upper_count: usize,
    lower_count: usize,
    balance: bool,
) -> (usize, usize, Option<String>) {
    let count_for = |area: f64, apl: f64| (area / apl).round().max(1.0) as usize;
    let (upper, lower) = (upper_count, lower_count);
    let upper_area = upper as f64 * upper_apl;
    let lower_area = lower as f64 * lower_apl;
    let mismatch = (upper_area - lower_area).abs() / upper_area.max(lower_area);
    if mismatch < 0.02 {
        return (upper, lower, None);
    }

    let heavy_area = upper_area.max(lower_area);
    let upper_light = upper_area < lower_area;
    if balance {
        if upper_light {
            let grown = count_for(heavy_area, upper_apl);
            let note = format!(
                "balanced: upper leaflet {upper} -> {grown} lipids to match the lower leaflet"
            );
            (grown, lower, Some(note))
        } else {
            let grown = count_for(heavy_area, lower_apl);
            let note = format!(
                "balanced: lower leaflet {lower} -> {grown} lipids to match the upper leaflet"
            );
            (upper, grown, Some(note))
        }
    } else {
        let (name, key, count_now, apl) = if upper_light {
            ("upper", "up", upper, upper_apl)
        } else {
            ("lower", "lo", lower, lower_apl)
        };
        let to_match = count_for(heavy_area, apl);
        let pct = (mismatch * 100.0).round();
        let note = format!(
            "leaflet areas differ by {pct}%: the {name} leaflet will start under-dense \
             (a possible pore). Set {key}={to_match} in --count (from {count_now}), \
             or drop --no-balance to grow it to match automatically"
        );
        (upper, lower, Some(note))
    }
}

/// Per-leaflet lipid composition. Each lipid molecule is classified into one of
/// the two leaflets and counted per species (in the system's species order): a
/// flat membrane splits by centroid z about the mean lipid bead z (`upper`/`lower`),
/// a vesicle by 3D centroid radius (`outer`/`inner`), and a tubule by 2D radial
/// distance in the x/y plane (`outer`/`inner`). Protein particles are excluded.
fn leaflet_split(sys: &EmbeddedSystem, geometry: Geometry) -> Vec<report::Leaflet> {
    let lipids = sys.lipid_particles();
    if lipids.is_empty() {
        return Vec::new();
    }
    let radial = matches!(geometry, Geometry::Vesicle | Geometry::Tubule);
    let n = lipids.len() as f64;
    let center = if radial {
        let mut c = [0.0; 3];
        for p in lipids {
            c[0] += p.position[0];
            c[1] += p.position[1];
            c[2] += p.position[2];
        }
        [c[0] / n, c[1] / n, c[2] / n]
    } else {
        [0.0; 3]
    };
    let scalar = |pos: [f64; 3]| match geometry {
        Geometry::Membrane => pos[2],
        Geometry::Vesicle => {
            let dx = pos[0] - center[0];
            let dy = pos[1] - center[1];
            let dz = pos[2] - center[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        }
        Geometry::Tubule => (pos[0] - center[0]).hypot(pos[1] - center[1]),
    };

    let threshold = lipids.iter().map(|p| scalar(p.position)).sum::<f64>() / n;
    let n_species = sys.lipid_counts.len();
    let mut far = vec![0usize; n_species];
    let mut near = vec![0usize; n_species];
    let mut species_idx = 0usize;
    for mol in lipids.chunk_by(|a, b| a.molecule_id == b.molecule_id) {
        let res = mol[0].residue.as_str();
        if species_idx >= n_species || sys.lipid_counts[species_idx].0 != res {
            if let Some(idx) = sys.lipid_counts.iter().position(|(s, _)| s == res) {
                species_idx = idx;
            } else {
                continue;
            }
        }
        let mean = mol.iter().map(|p| scalar(p.position)).sum::<f64>() / mol.len() as f64;
        if mean >= threshold {
            far[species_idx] += 1;
        } else {
            near[species_idx] += 1;
        }
    }

    let in_species_order = |counts: &[usize]| -> Vec<(String, usize)> {
        sys.lipid_counts
            .iter()
            .zip(counts)
            .filter(|&(_, &c)| c > 0)
            .map(|((species, _), &c)| (species.clone(), c))
            .collect()
    };
    let (far_label, near_label) = if radial {
        ("outer", "inner")
    } else {
        ("upper", "lower")
    };
    vec![
        (far_label, in_species_order(&far)),
        (near_label, in_species_order(&near)),
    ]
}

#[allow(clippy::too_many_arguments)]
fn write_outputs(
    system: &mut EmbeddedSystem,
    force_field: &ForceField,
    common: &CommonBuildArgs,
    geometry: Geometry,
    center: CenterArg,
    controls: &ControlSpecs,
    base: &str,
    gif_frames: &[RelaxFrame],
    gif_bounds: [f64; 4],
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&common.out)?;
    let named = |extension: &str| common.out.join(format!("{base}.{extension}"));

    system.center_in_box();
    let box_size = system.box_size;
    if matches!(center, CenterArg::Origin) {
        system.translate([-box_size[0] / 2.0, -box_size[1] / 2.0, -box_size[2] / 2.0]);
    }
    let system: &EmbeddedSystem = system;

    let want_inp = common.want_inp();
    let want_vmd = common.want_vmd();
    let want_gif = common.want_gif();

    let formats: [StructureFormat; 5] = [
        (true, "gro", EmbeddedSystem::to_gro),
        (common.pdb, "pdb", EmbeddedSystem::to_pdb),
        (common.psf || want_vmd, "psf", EmbeddedSystem::to_psf),
        (common.crd, "crd", EmbeddedSystem::to_crd),
        (common.cif, "cif", EmbeddedSystem::to_cif),
    ];
    for (enabled, ext, render) in formats {
        if enabled {
            write(&named(ext), &render(system))?;
        }
    }

    if common.top || want_inp {
        write(&named("top"), &system.to_top())?;
        write(&common.out.join("isolf.itp"), &force_field.to_itp())?;
        if system.has_proteins() {
            write(&common.out.join("aicg2p.itp"), &render_aicg2p_itp())?;
            for tmpl in &system.protein_templates {
                write(
                    &common.out.join(format!("{}.itp", tmpl.name)),
                    &tmpl.to_itp(),
                )?;
            }
        }
    }

    if want_gif && !gif_frames.is_empty() {
        let opts = anim::GifOptions {
            mode: match common.gif_mode.unwrap_or(GifModeArg::Point) {
                GifModeArg::Point => anim::GifMode::Point,
                GifModeArg::Density => anim::GifMode::Density,
            },
            scale: common.gif_scale.unwrap_or(1.0),
            fps: common.gif_fps.unwrap_or(14),
        };
        let path = named("gif");
        anim::write_gif(&path, gif_frames, gif_bounds, &opts)?;
        report::record(&path);
    }

    if want_inp {
        write_control_files(&common.out, base, geometry, box_size, controls, rng)?;
    }
    if want_vmd {
        write_vmd_scripts(
            &common.out,
            base,
            geometry,
            center,
            system,
            force_field,
            common,
            controls,
        )?;
    }
    Ok(())
}

/// Write the GENESIS control files: minimize, then equilibrate and run
/// production under the build mode's ensemble (NPT for a flat membrane, NVT for
/// a vesicle).
fn write_control_files(
    out: &Path,
    base: &str,
    geometry: Geometry,
    box_size: [f64; 3],
    controls: &ControlSpecs,
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    let min_params = controls.min.to_params(box_size, controls.temperature);
    write(
        &out.join("min.inp"),
        &isolf::inp::minimization(base, &min_params)?,
    )?;

    let prod_seed = rng.random_range(1..=MAX_ISEED);
    let equil_seed = rng.random_range(1..=MAX_ISEED);
    let production = controls.prod.to_params(controls.temperature, prod_seed);
    let equilibration = controls.equil.to_params(controls.temperature, equil_seed);

    if geometry.is_npt() {
        write(
            &out.join("npt.inp"),
            &isolf::inp::npt(base, "min.rst", "npt", &equilibration)?,
        )?;
        write(
            &out.join("pro.inp"),
            &isolf::inp::npt(base, "npt.rst", "pro", &production)?,
        )?;
    } else {
        write(
            &out.join("nvt.inp"),
            &isolf::inp::nvt(base, "min.rst", "nvt", &equilibration)?,
        )?;
        write(
            &out.join("pro.inp"),
            &isolf::inp::nvt(base, "nvt.rst", "pro", &production)?,
        )?;
    }
    Ok(())
}

/// Write VMD visualization scripts for the built structure (and per MD stage
/// when `--inp` is enabled).
#[allow(clippy::too_many_arguments)]
fn write_vmd_scripts(
    out: &Path,
    base: &str,
    geometry: Geometry,
    center: CenterArg,
    system: &EmbeddedSystem,
    force_field: &ForceField,
    common: &CommonBuildArgs,
    controls: &ControlSpecs,
) -> Result<(), Box<dyn Error>> {
    let background = common.vmd_background.as_deref().unwrap_or("white");
    let material = common.vmd_material.as_deref().unwrap_or("AOChalky");
    let options = VmdOptions {
        surface: common.vmd_surface,
        cutaway: common.vmd_cutaway,
        render: common.vmd_render,
        background,
        material,
        smooth: common.vmd_smooth.unwrap_or(5),
        lipid_coloring: match common.vmd_color_by.unwrap_or(LipidColorScope::Head) {
            LipidColorScope::Head => isolf::vmd::LipidColoring::Head,
            LipidColorScope::Bead => isolf::vmd::LipidColoring::Role,
            LipidColorScope::Whole => isolf::vmd::LipidColoring::Whole,
        },
        recenter_to_box: matches!(center, CenterArg::Origin),
    };
    let layout = geometry.layout();
    let script = |source| {
        isolf::vmd::vmd_script(
            force_field,
            &system.particles,
            base,
            layout,
            source,
            &options,
        )
    };

    write(
        &out.join(format!("{base}.vmd")),
        &script(VmdSource::Structure),
    )?;

    if common.want_inp() {
        let equil_name = if geometry.is_npt() { "npt" } else { "nvt" };
        let stages = [
            ("min", controls.min.num_steps, controls.min.crdout_period),
            (
                equil_name,
                controls.equil.num_steps,
                controls.equil.crdout_period,
            ),
            ("pro", controls.prod.num_steps, controls.prod.crdout_period),
        ];
        for (name, num_steps, crdout_period) in stages {
            let frames = (crdout_period > 0).then(|| num_steps / crdout_period);
            write(
                &out.join(format!("{name}.vmd")),
                &script(VmdSource::Trajectory { name, frames }),
            )?;
        }
    }
    Ok(())
}

fn write(path: &Path, contents: &str) -> std::io::Result<()> {
    fs::write(path, contents)?;
    report::record(path);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn balanced_counts_covers_the_leaflet_cases() {
        let (u, l, note) = balanced_counts(0.64, 0.64, 200, 200, false);
        assert_eq!((u, l), (200, 200));
        assert!(note.is_none());

        let (u, l, note) = balanced_counts(0.64, 0.64, 200, 150, false);
        assert_eq!((u, l), (200, 150));
        assert!(note.unwrap().contains("lo=200 in --count"));

        let (u, l, note) = balanced_counts(0.64, 0.64, 200, 150, true);
        assert_eq!((u, l), (200, 200));
        assert!(note.unwrap().contains("balanced"));
    }

    #[test]
    fn parses_geometry_specs() {
        assert_eq!(parse_count_spec(&["256".into()]).unwrap(), Some((256, 256)));
        assert_eq!(
            parse_count_spec(&["up=300,lo=200".into()]).unwrap(),
            Some((300, 200))
        );
        assert!(parse_count_spec(&["foo=10".into()]).is_err());

        assert_eq!(parse_box_spec(&["15".into()]).unwrap(), Some((15.0, 15.0)));
        assert_eq!(
            parse_box_spec(&["x=10".into(), "y=20".into()]).unwrap(),
            Some((10.0, 20.0))
        );
        assert!(parse_box_spec(&["z=10".into()]).is_err());

        assert_eq!(
            parse_radius_spec(&["20".into()]).unwrap(),
            Some(VesicleRadius::Outer(20.0))
        );
        assert_eq!(
            parse_radius_spec(&["ri=14,ro=20".into()]).unwrap(),
            Some(VesicleRadius::Both {
                inner: 14.0,
                outer: 20.0
            })
        );
        assert!(parse_radius_spec(&["r=20".into()]).is_err());
    }

    #[test]
    fn parses_stage_override_specs() {
        let m = parse_min_spec(&["2000".into(), "eneout=20,crdout=500".into()]).unwrap();
        assert_eq!(m.num_steps, 2000);
        assert_eq!(m.eneout_period, 20);
        assert_eq!(m.crdout_period, 500);
        assert_eq!(m.nbupdate_period, 10);
        assert!(parse_min_spec(&["dt=0.01".into()]).is_err());

        let d = parse_dyn_spec(
            "--equil",
            &["steps=50000,dt=0.015,temp=310,seed=7".into()],
            DynSpec::default_equil(),
        )
        .unwrap();
        assert_eq!(d.num_steps, 50_000);
        assert_eq!(d.time_step, 0.015);
        assert_eq!(d.temperature, Some(310.0));
        assert_eq!(d.seed, Some(7));
        assert!(parse_dyn_spec("--prod", &["bad=1".into()], DynSpec::default_prod()).is_err());
    }

    #[test]
    fn suggests_close_lipid_names() {
        let ff = ForceField::isolf();
        let hits = suggest_lipids(&ff, "POPX");
        assert!(hits.contains(&"POPC"));
        assert!(hits.contains(&"POPS"));
        assert!(suggest_lipids(&ff, "ZZZZZZ").is_empty());
    }

    #[test]
    fn cli_subcommands_and_auto_enabled_flags() {
        // Symmetric membrane with auto-enabled --inp, --vmd, --gif.
        let cli = Cli::try_parse_from([
            "isolf",
            "membrane",
            "-L",
            "POPC=3,POPS=1",
            "-n",
            "256",
            "--equil",
            "steps=50000",
            "--vmd-cutaway",
            "--gif-mode",
            "density",
        ])
        .unwrap();
        let Command::Membrane(args) = cli.command else {
            panic!("expected membrane subcommand");
        };
        assert!(args.common.want_inp());
        assert!(args.common.want_vmd());
        assert!(args.common.want_gif());

        // --lipids conflicts with --upper.
        assert!(
            Cli::try_parse_from([
                "isolf", "membrane", "-L", "POPC", "--upper", "DOPC", "-n", "64"
            ])
            .is_err()
        );

        // Membrane sized by --box (including combined --count + --box).
        let cli = Cli::try_parse_from([
            "isolf",
            "membrane",
            "-L",
            "POPC",
            "--box",
            "x=20,y=10",
            "--count",
            "100",
        ])
        .unwrap();
        assert!(matches!(cli.command, Command::Membrane(_)));

        // Vesicle with asymmetric shells.
        let cli = Cli::try_parse_from([
            "isolf",
            "vesicle",
            "--outer",
            "POPC",
            "--inner",
            "POPS",
            "-r",
            "ri=14,ro=20",
        ])
        .unwrap();
        assert!(matches!(cli.command, Command::Vesicle(_)));

        // Tubule with phase-separated raft and anti-registration.
        let cli = Cli::try_parse_from([
            "isolf",
            "tubule",
            "-L",
            "DPPC=1,DSPC=1;POPC=2,POPS=1",
            "-r",
            "10",
            "-l",
            "25",
            "--phase",
            "raft",
            "--anti-register",
        ])
        .unwrap();
        let Command::Tubule(tub) = cli.command else {
            panic!("expected tubule subcommand");
        };
        assert_eq!(tub.length, 25.0);
        assert_eq!(tub.common.phase, PhaseArg::Raft);
        assert!(tub.common.anti_register);

        // Standalone protein subcommand with IDR, Clementi option, and auto-enabled --inp/--vmd.
        let cli = Cli::try_parse_from([
            "isolf",
            "protein",
            "channel.pdb",
            "--model",
            "clementi",
            "--idr",
            "1-20",
            "--vmd-surface",
            "--equil",
            "steps=10000",
        ])
        .unwrap();
        let Command::Protein(prot) = cli.command else {
            panic!("expected protein subcommand");
        };
        assert_eq!(prot.model, ProteinModelArg::ClementiGo);
        assert!(prot.want_vmd());
        assert!(prot.want_inp());

        // Membrane with repeatable -p protein specs.
        let cli = Cli::try_parse_from([
            "isolf",
            "membrane",
            "-L",
            "POPC",
            "-b",
            "20",
            "-p",
            "a.pdb,z=center,copies=2",
            "-p",
            "b.cif,x=3,tilt=15",
        ])
        .unwrap();
        let Command::Membrane(mem) = cli.command else {
            panic!("expected membrane subcommand");
        };
        assert_eq!(mem.common.proteins.len(), 2);

        // Discovery subcommand.
        let cli = Cli::try_parse_from(["isolf", "lipids"]).unwrap();
        assert!(matches!(cli.command, Command::Lipids(_)));
    }
}
