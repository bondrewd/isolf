//! Clap argument definitions and enum adapters for every `isolf` subcommand.

use std::error::Error;
use std::path::PathBuf;

use clap::{Parser, Subcommand};

use isolf::inp::Ensemble;
use isolf::membrane::{DEFAULT_PADDING, DEFAULT_TEMPERATURE, Phase};
use isolf::protein::aicg2p::{AicgScale, ProteinModel, SafeDihedral};
use isolf::umbrella::{DEFAULT_UMBRELLA_BOX_Z, DEFAULT_UMBRELLA_K};

use super::parse::{ControlSpecs, DynSpec, parse_dyn_spec, parse_min_spec};

/// Where the finished system is centred (`--center`). The periodic box is the
/// same either way; only the coordinates move.
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum CenterArg {
    /// Centred on the coordinate origin (0, 0, 0) (default).
    #[default]
    Origin,
    /// Centred in the periodic box [0, box].
    Box,
}

/// Scope for `--vmd-color-by`: how lipid beads are coloured.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LipidColorScope {
    /// Colour the head bead by lipid species; the rest by bead role.
    Head,
    /// Colour every bead by its role (head/phosphate/glycerol/tail).
    Bead,
    /// Colour every bead by lipid species.
    Whole,
}

/// Scope for `--gif-mode`: how `--gif` renders each frame.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum GifModeArg {
    /// One disc per lipid; clear for small to medium systems.
    Point,
    /// Local-density heatmap; better for large, dense systems.
    Density,
}

/// Initial lateral lipid distribution (`--phase`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum PhaseArg {
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
pub(crate) enum ProteinCenterArg {
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
pub(crate) enum ProteinModelArg {
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
pub(crate) enum AicgScaleArg {
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
pub(crate) enum SafeDihedralArg {
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

/// Umbrella sampling execution mode (`isolf umbrella setup --mode`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum UmbrellaModeArg {
    /// Independent per-window directories (win_1/, win_2/, ...) (default).
    #[default]
    Windows,
    /// Single multi-replica GENESIS REUS ([REMD] type1 = RESTRAINT) setup.
    Reus,
}

/// Dynamics ensemble for umbrella sampling (`isolf umbrella setup --ensemble`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum UmbrellaEnsembleArg {
    /// Constant volume NVT (recommended for fixed-area bilayer PMFs, default).
    #[default]
    Nvt,
    /// Constant pressure NPT (Z-FIXED-SEMI-ISO).
    Npt,
}

impl UmbrellaEnsembleArg {
    pub(crate) fn to_ensemble(self) -> Ensemble {
        match self {
            Self::Nvt => Ensemble::Nvt,
            Self::Npt => Ensemble::Npt,
        }
    }

    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Nvt => "nvt",
            Self::Npt => "npt",
        }
    }
}

/// GENESIS restraint function (`isolf umbrella setup --restraint-func`).
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum RestraintFuncArg {
    /// Center-of-mass distance (`DISTMASS`, default).
    #[default]
    #[value(name = "distmass")]
    DistMass,
    /// Geometric-center distance (`DIST`).
    #[value(name = "dist")]
    Dist,
}

impl RestraintFuncArg {
    pub(crate) fn as_genesis_str(self) -> &'static str {
        match self {
            Self::DistMass => "DISTMASS",
            Self::Dist => "DIST",
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
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) command: Command,
}

#[derive(Subcommand, Debug)]
pub(crate) enum Command {
    /// Build a flat lipid bilayer.
    Membrane(MembraneArgs),
    /// Build a spherical lipid vesicle.
    Vesicle(VesicleArgs),
    /// Build a cylindrical lipid tubule (periodic along z).
    Tubule(TubuleArgs),
    /// Coarse-grain a protein PDB/mmCIF with AICG2+ or Clementi Gō for GENESIS.
    Protein(ProteinArgs),
    /// Set up 1D bilayer umbrella sampling (US / REUS) and extract CVs for MBAR.
    #[command(subcommand)]
    Umbrella(UmbrellaCommand),
    /// Generate or update lipid_protein.itp from parameter files or inline overrides.
    Params(ParamsArgs),
    /// List all supported iSoLF lipids.
    Lipids(LipidsArgs),
    /// Update isolf in place to the latest release.
    Update,
    /// Uninstall isolf, removing the binary from the system.
    Uninstall,
}

/// Subcommands under `isolf umbrella`.
#[derive(Subcommand, Debug)]
pub(crate) enum UmbrellaCommand {
    /// Build umbrella-sampling windows (or REUS) + mbar.inp for amino acids or PDB/mmCIF proteins.
    Setup(Box<UmbrellaSetupArgs>),
    /// Extract RESTR_CVS001 from GENESIS atdyn logs (and demux REUS .rem files) into cv/{}.dat.
    Cv(UmbrellaCvArgs),
}

/// Flat-bilayer composition and sizing flags shared between `membrane` and `umbrella setup`.
#[derive(clap::Args, Debug)]
pub(crate) struct BilayerArgs {
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
        help_heading = "Bilayer"
    )]
    pub(crate) lipids: Vec<String>,

    /// Upper leaflet composition as NAME=WEIGHT pairs (mirrored to lower if --lower is omitted).
    #[arg(short = 'u', long, value_name = "SPEC", help_heading = "Bilayer")]
    pub(crate) upper: Vec<String>,

    /// Lower leaflet composition as NAME=WEIGHT pairs (mirrored to upper if --upper is omitted).
    #[arg(short = 'l', long, value_name = "SPEC", help_heading = "Bilayer")]
    pub(crate) lower: Vec<String>,

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
        help_heading = "Bilayer"
    )]
    pub(crate) count: Vec<String>,

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
        help_heading = "Bilayer"
    )]
    pub(crate) box_size: Vec<String>,
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
pub(crate) struct MembraneArgs {
    #[command(flatten)]
    pub(crate) bilayer: BilayerArgs,

    /// System centre: origin (0,0,0) or box (Lx/2,Ly/2,Lz/2).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = CenterArg::Origin,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    pub(crate) center: CenterArg,

    /// Vacuum padding in nm added on each side along z.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    pub(crate) padding: f64,

    /// Keep unequal leaflet counts as given instead of growing the lighter leaflet to match area.
    #[arg(long, help_heading = "Geometry")]
    pub(crate) no_balance: bool,

    #[command(flatten)]
    pub(crate) common: CommonBuildArgs,
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
pub(crate) struct VesicleArgs {
    /// Lipid composition for both shells (symmetric vesicle) as NAME=WEIGHT pairs.
    #[arg(
        short = 'L',
        long,
        value_name = "SPEC",
        conflicts_with_all = ["outer", "inner"],
        help_heading = "Lipid composition"
    )]
    pub(crate) lipids: Vec<String>,

    /// Outer shell composition as NAME=WEIGHT pairs (copied to inner if --inner is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    pub(crate) outer: Vec<String>,

    /// Inner (lumen) shell composition as NAME=WEIGHT pairs (copied to outer if --outer is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    pub(crate) inner: Vec<String>,

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
    pub(crate) radius: Vec<String>,

    /// System centre: origin (0,0,0) or box (Lx/2,Ly/2,Lz/2).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = CenterArg::Origin,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    pub(crate) center: CenterArg,

    /// Vacuum padding in nm added on every side of the cubic box.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    pub(crate) padding: f64,

    #[command(flatten)]
    pub(crate) common: CommonBuildArgs,
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
pub(crate) struct TubuleArgs {
    /// Lipid composition for both shells (symmetric tubule) as NAME=WEIGHT pairs.
    #[arg(
        short = 'L',
        long,
        value_name = "SPEC",
        conflicts_with_all = ["outer", "inner"],
        help_heading = "Lipid composition"
    )]
    pub(crate) lipids: Vec<String>,

    /// Outer shell composition as NAME=WEIGHT pairs (copied to inner if --inner is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    pub(crate) outer: Vec<String>,

    /// Inner (lumen) shell composition as NAME=WEIGHT pairs (copied to outer if --outer is omitted).
    #[arg(long, value_name = "SPEC", help_heading = "Lipid composition")]
    pub(crate) inner: Vec<String>,

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
    pub(crate) radius: Vec<String>,

    /// Periodic axial length of the tubule along z in nm.
    #[arg(
        short = 'l',
        long,
        required = true,
        value_name = "NM",
        help_heading = "Geometry"
    )]
    pub(crate) length: f64,

    /// System centre: origin (0,0,0) or box (Lx/2,Ly/2,Lz/2).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = CenterArg::Origin,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    pub(crate) center: CenterArg,

    /// Vacuum padding in nm added on each side along x and y.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    pub(crate) padding: f64,

    #[command(flatten)]
    pub(crate) common: CommonBuildArgs,
}

/// Options shared between `membrane`, `vesicle`, and `tubule`.
#[derive(clap::Args, Debug)]
pub(crate) struct CommonBuildArgs {
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
    pub(crate) phase: PhaseArg,

    /// Offset phase domains on opposite leaflets/shells (anti-registered domains).
    #[arg(long, help_heading = "Lipid composition")]
    pub(crate) anti_register: bool,

    /// Embed a PDB/mmCIF protein or single AA bead (repeatable) as SPEC[,KEY=VAL,...].
    ///
    /// Preserves input orientation by default (OPM midplane z = 0).
    /// Keys: z=<nm|center>, x=<nm>, y=<nm>, tilt=<deg>, rot=<deg>, copies=<N>,
    ///       clearance=<nm>, charge=<e>, model=<aicg2+|clementi>, scale=<general|average>,
    ///       safe-dihedral=<0|1|3>, idr=<ranges>, idr-local=<ranges>,
    ///       idr-nonlocal=<ranges>, charges=<file>, name=<molname>
    /// Examples:
    ///   -p channel.pdb,z=center,copies=2,clearance=0.2
    ///   -p aa:TRP,z=1.5
    #[arg(
        short = 'p',
        long = "protein",
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Proteins"
    )]
    pub(crate) proteins: Vec<String>,

    /// Lipid–protein interaction parameter file; writes lipid_protein.itp and includes it in .top.
    #[arg(long, value_name = "FILE", help_heading = "Proteins")]
    pub(crate) lipid_protein: Option<PathBuf>,

    /// Output directory.
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output files"
    )]
    pub(crate) out: PathBuf,

    /// Base name for output files [default: membrane / vesicle / tubule].
    #[arg(long, value_name = "NAME", help_heading = "Output files")]
    pub(crate) name: Option<String>,

    /// Write topology (.top) and force-field (isolf.itp) files.
    #[arg(long, help_heading = "Output files")]
    pub(crate) top: bool,

    /// Write GENESIS .inp control files for min, equil, and prod (implies --top).
    #[arg(long, help_heading = "Output files")]
    pub(crate) inp: bool,

    /// Write VMD .vmd visualization scripts (implies --psf).
    #[arg(long, help_heading = "Output files")]
    pub(crate) vmd: bool,

    /// Write a .gif animating the leaflet relaxation.
    #[arg(long, help_heading = "Output files")]
    pub(crate) gif: bool,

    #[command(flatten)]
    pub(crate) formats: FormatArgs,

    #[command(flatten)]
    pub(crate) controls: ControlArgs,

    #[command(flatten)]
    pub(crate) vmd_opts: VmdArgs,

    /// Colouring scheme [default: head] (implies --vmd).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        value_name = "SCOPE",
        help_heading = "VMD visualization (--vmd)"
    )]
    pub(crate) vmd_color_by: Option<LipidColorScope>,

    /// Render mode: point scatter or density heatmap [default: point] (implies --gif).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        value_name = "MODE",
        help_heading = "Relaxation animation (--gif)"
    )]
    pub(crate) gif_mode: Option<GifModeArg>,

    /// Resolution multiplier, e.g. 2 doubles it [default: 1] (implies --gif).
    #[arg(
        long,
        value_name = "FACTOR",
        help_heading = "Relaxation animation (--gif)"
    )]
    pub(crate) gif_scale: Option<f64>,

    /// Playback speed in frames per second [default: 14] (implies --gif).
    #[arg(
        long,
        value_name = "FPS",
        help_heading = "Relaxation animation (--gif)"
    )]
    pub(crate) gif_fps: Option<u16>,

    #[command(flatten)]
    pub(crate) general: GeneralArgs,
}

impl CommonBuildArgs {
    pub(crate) fn want_inp(&self) -> bool {
        self.inp || self.controls.has_overrides()
    }

    pub(crate) fn want_vmd(&self) -> bool {
        self.vmd || self.vmd_opts.any_set() || self.vmd_color_by.is_some()
    }

    pub(crate) fn want_gif(&self) -> bool {
        self.gif || self.gif_mode.is_some() || self.gif_scale.is_some() || self.gif_fps.is_some()
    }
}

/// Optional coordinate / structure formats shared across all build subcommands.
#[derive(clap::Args, Debug)]
pub(crate) struct FormatArgs {
    /// Write a .pdb coordinate file.
    #[arg(long, help_heading = "Output files")]
    pub(crate) pdb: bool,

    /// Write a CHARMM .psf structure file.
    #[arg(long, help_heading = "Output files")]
    pub(crate) psf: bool,

    /// Write a CHARMM .crd coordinate file.
    #[arg(long, help_heading = "Output files")]
    pub(crate) crd: bool,

    /// Write an mmCIF .cif coordinate file.
    #[arg(long, help_heading = "Output files")]
    pub(crate) cif: bool,
}

/// GENESIS `.inp` stage override flags shared across all build subcommands.
#[derive(clap::Args, Debug)]
pub(crate) struct ControlArgs {
    /// Simulation temperature in K for GENESIS .inp files [default: 323.15].
    #[arg(
        short = 't',
        long,
        value_name = "K",
        help_heading = "GENESIS control files (--inp)"
    )]
    pub(crate) temperature: Option<f64>,

    /// Explicit cg_cutoffdist_ele in Angstrom (pairlist = cutoff + 5 A); auto-scaled for small boxes if unset.
    #[arg(long, value_name = "A", help_heading = "GENESIS control files (--inp)")]
    pub(crate) ele_cutoff: Option<f64>,

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
    pub(crate) min: Vec<String>,

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
    pub(crate) equil: Vec<String>,

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
    pub(crate) prod: Vec<String>,
}

impl ControlArgs {
    pub(crate) fn has_overrides(&self) -> bool {
        !self.min.is_empty() || !self.equil.is_empty() || !self.prod.is_empty()
    }

    /// Parse and pre-validate `--min`, `--equil`, `--prod`, `--temperature`,
    /// and `--ele-cutoff` before building so bad period multiples fail immediately.
    pub(crate) fn prepare(
        &self,
        want_inp: bool,
        is_npt: bool,
    ) -> Result<(ControlSpecs, Option<String>), Box<dyn Error>> {
        if let Some(t) = self.temperature
            && (!t.is_finite() || t <= 0.0)
        {
            return Err(format!("--temperature: '{t}' is not a positive temperature in K").into());
        }
        if let Some(rc) = self.ele_cutoff
            && (!rc.is_finite() || rc <= 0.0)
        {
            return Err(format!("--ele-cutoff: '{rc}' must be a positive distance in A").into());
        }
        let temperature = self.temperature.unwrap_or(DEFAULT_TEMPERATURE);
        let min = parse_min_spec(&self.min)?;
        let equil = parse_dyn_spec("--equil", &self.equil, DynSpec::default_equil())?;
        let prod = parse_dyn_spec("--prod", &self.prod, DynSpec::default_prod())?;

        if want_inp {
            isolf::inp::minimization("x", &min.to_params([20.0; 3], temperature, self.ele_cutoff))
                .map_err(|e| format!("--min: {e}"))?;
            let check_dyn = |flag: &str, spec: DynSpec| -> Result<(), Box<dyn Error>> {
                let p = spec.to_params(temperature, 1, [20.0; 3], self.ele_cutoff);
                let res = if is_npt {
                    isolf::inp::npt("x", "x", "x", &p)
                } else {
                    isolf::inp::nvt("x", "x", "x", &p)
                };
                res.map(|_| ()).map_err(|e| format!("{flag}: {e}").into())
            };
            check_dyn("--equil", equil)?;
            check_dyn("--prod", prod)?;
        }

        let temp_note = (self.temperature.is_some() && !want_inp).then(|| {
            "--temperature only affects GENESIS .inp files (pass --inp to write them)".to_string()
        });
        Ok((
            ControlSpecs {
                temperature,
                ele_cutoff: self.ele_cutoff,
                min,
                equil,
                prod,
            },
            temp_note,
        ))
    }
}

/// VMD visualization options shared across all build subcommands.
#[derive(clap::Args, Debug)]
pub(crate) struct VmdArgs {
    /// Use a smooth QuickSurf surface instead of VDW spheres / CA Tube (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    pub(crate) vmd_surface: bool,

    /// Slice the system with a clipping plane to reveal the cross-section (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    pub(crate) vmd_cutaway: bool,

    /// Append a Tachyon ray-tracing render command to each script (implies --vmd).
    #[arg(long, help_heading = "VMD visualization (--vmd)")]
    pub(crate) vmd_render: bool,

    /// VMD background colour name [default: white] (implies --vmd).
    #[arg(long, value_name = "COLOR", help_heading = "VMD visualization (--vmd)")]
    pub(crate) vmd_background: Option<String>,

    /// VMD representation material, e.g. AOChalky or Goodsell [default: AOChalky] (implies --vmd).
    #[arg(long, value_name = "NAME", help_heading = "VMD visualization (--vmd)")]
    pub(crate) vmd_material: Option<String>,

    /// Trajectory smoothing window in frames [default: 5] (implies --vmd).
    #[arg(long, value_name = "N", help_heading = "VMD visualization (--vmd)")]
    pub(crate) vmd_smooth: Option<usize>,
}

impl VmdArgs {
    pub(crate) fn any_set(&self) -> bool {
        self.vmd_surface
            || self.vmd_cutaway
            || self.vmd_render
            || self.vmd_background.is_some()
            || self.vmd_material.is_some()
            || self.vmd_smooth.is_some()
    }
}

/// Seed and terminal-output flags shared across all build subcommands.
#[derive(clap::Args, Debug)]
pub(crate) struct GeneralArgs {
    /// Seed value for the RNG (random if unset).
    #[arg(long, value_name = "SEED", help_heading = "General")]
    pub(crate) seed: Option<u64>,

    /// Print only the final summary and any warnings.
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        help_heading = "General"
    )]
    pub(crate) quiet: bool,

    /// Print full file paths in the Write section.
    #[arg(short = 'v', long, help_heading = "General")]
    pub(crate) verbose: bool,

    /// Disable coloured terminal output.
    #[arg(long, help_heading = "General")]
    pub(crate) no_color: bool,

    /// Use ASCII status markers instead of Unicode glyphs.
    #[arg(long, help_heading = "General")]
    pub(crate) ascii: bool,
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
pub(crate) struct ProteinArgs {
    /// Input all-atom protein structure (.pdb or .cif/.mmcif).
    #[arg(value_name = "INPUT")]
    pub(crate) input: PathBuf,

    /// Coarse-grained protein force field: aicg2+ (default) or clementi.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = ProteinModelArg::Aicg2Plus,
        value_name = "MODEL",
        help_heading = "Protein model"
    )]
    pub(crate) model: ProteinModelArg,

    /// AICG2+ energy scaling scheme: general (default) or average.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = AicgScaleArg::General,
        value_name = "SCHEME",
        help_heading = "Protein model"
    )]
    pub(crate) scale: AicgScaleArg,

    /// Safe dihedral potential mode: 0 (none), 1 (cos2, default), or 3 (sin3).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = SafeDihedralArg::Cos2,
        value_name = "MODE",
        help_heading = "Protein model"
    )]
    pub(crate) safe_dihedral: SafeDihedralArg,

    /// Disordered residue ranges excluded from both local and nonlocal Gō terms (e.g. "1-25,80-95").
    #[arg(long, value_name = "RANGES", help_heading = "Protein model")]
    pub(crate) idr: Option<String>,

    /// Disordered residue ranges excluded from AICG2+ local Gaussian 1-3 and 1-4 terms.
    #[arg(long, value_name = "RANGES", help_heading = "Protein model")]
    pub(crate) idr_local: Option<String>,

    /// Disordered residue ranges excluded from native Gō contacts ([ pairs ]).
    #[arg(long, value_name = "RANGES", help_heading = "Protein model")]
    pub(crate) idr_nonlocal: Option<String>,

    /// RESPAC / user charge override table file (<1-based index> <charge> per line).
    #[arg(long, value_name = "FILE", help_heading = "Protein model")]
    pub(crate) charges: Option<PathBuf>,

    /// Native contact energy scale factor for the Clementi Gō model.
    #[arg(
        long,
        value_name = "FACTOR",
        default_value_t = 1.0,
        help_heading = "Protein model"
    )]
    pub(crate) contact_scale: f64,

    /// Coordinate centering: keep (input coordinates, default), origin, or box.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = ProteinCenterArg::Keep,
        value_name = "WHERE",
        help_heading = "Geometry"
    )]
    pub(crate) center: ProteinCenterArg,

    /// Vacuum padding in nm added on each side of the periodic box.
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_PADDING,
        help_heading = "Geometry"
    )]
    pub(crate) padding: f64,

    /// Output directory.
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output files"
    )]
    pub(crate) out: PathBuf,

    /// Base / moleculetype name for output files [default: <stem>_cg].
    #[arg(long, value_name = "NAME", help_heading = "Output files")]
    pub(crate) name: Option<String>,

    /// Write GENESIS .inp control files for min, nvt, and pro.
    #[arg(long, help_heading = "Output files")]
    pub(crate) inp: bool,

    /// Write VMD .vmd visualization scripts (implies --psf).
    #[arg(long, help_heading = "Output files")]
    pub(crate) vmd: bool,

    #[command(flatten)]
    pub(crate) formats: FormatArgs,

    #[command(flatten)]
    pub(crate) controls: ControlArgs,

    #[command(flatten)]
    pub(crate) vmd_opts: VmdArgs,

    #[command(flatten)]
    pub(crate) general: GeneralArgs,
}

impl ProteinArgs {
    pub(crate) fn want_inp(&self) -> bool {
        self.inp || self.controls.has_overrides()
    }

    pub(crate) fn want_vmd(&self) -> bool {
        self.vmd || self.vmd_opts.any_set()
    }
}

/// Arguments for `isolf umbrella setup`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Build 1D bilayer umbrella-sampling windows (or REUS) and mbar.inp for GENESIS",
    long_about = "\
Place a single amino-acid bead (e.g. --probe TRP or --probe all) or a peptide /
protein PDB/mmCIF (e.g. --probe walp23.pdb,tilt=90) at z windows across a lipid
bilayer centered at z = 0 without carving lipids (preserving an identical
Hamiltonian across windows for MBAR).

Writes decoupled lipid_protein.itp ([ cg_ISOLF_parameters ] + head WCA/LJ),
GENESIS atdyn controls (DISTMASS + direction1 = Z2, U = k * (dz - z0)^2),
mbar.inp for mbar_analysis, and manifest.json."
)]
pub(crate) struct UmbrellaSetupArgs {
    /// Solute probe: AA code(s) ('TRP', 'TRP,ARG', 'all', 'aa:TRP,charge=0') or PDB/mmCIF[,KEY=VAL,...].
    ///
    /// PDB/mmCIF keys match `-p`: tilt=<deg>, rot=<deg>, model=<aicg2+|clementi>,
    /// scale=<general|average>, safe-dihedral=<0|1|3>, idr=<ranges>, charges=<file>, name=<molname>.
    #[arg(
        long,
        required = true,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Probe & parameters"
    )]
    pub(crate) probe: String,

    /// Optional charge override (e) for single-bead amino acid probe(s).
    #[arg(long, value_name = "Q", help_heading = "Probe & parameters")]
    pub(crate) charge: Option<f64>,

    /// Lipid–protein parameter file ([cg_ISOLF_parameters], optional [cg_WCA_parameters] / [cg_LJ_parameters]).
    #[arg(long, value_name = "FILE", help_heading = "Probe & parameters")]
    pub(crate) lipid_protein: Option<PathBuf>,

    /// Override probe [ cg_ISOLF_parameters ] epsilon in kJ/mol [default: 1.5].
    #[arg(long, value_name = "KJ/MOL", help_heading = "Probe & parameters")]
    pub(crate) eps: Option<f64>,

    /// Override probe [ cg_ISOLF_parameters ] sigma in nm [default: rmin / 2^(1/6)].
    #[arg(long, value_name = "NM", help_heading = "Probe & parameters")]
    pub(crate) sigma: Option<f64>,

    /// Override probe [ cg_ISOLF_parameters ] omega in nm [default: 1.0].
    #[arg(long, value_name = "NM", help_heading = "Probe & parameters")]
    pub(crate) omega: Option<f64>,

    /// Do not auto-generate Lorentz–Berthelot WCA repulsion between probe residues and lipid head beads.
    #[arg(long, help_heading = "Probe & parameters")]
    pub(crate) no_auto_wca: bool,

    /// Reuse a pre-equilibrated iSoLF bilayer .gro file instead of building a new bilayer.
    #[arg(
        long = "bilayer",
        value_name = "GRO",
        conflicts_with_all = ["lipids", "upper", "lower", "count", "box_size"],
        help_heading = "Bilayer"
    )]
    pub(crate) bilayer_gro: Option<PathBuf>,

    #[command(flatten)]
    pub(crate) bilayer: BilayerArgs,

    /// Box length along z in nm (expanded automatically if windows require more clearance).
    #[arg(
        long,
        value_name = "NM",
        default_value_t = DEFAULT_UMBRELLA_BOX_Z,
        help_heading = "Bilayer"
    )]
    pub(crate) box_z: f64,

    /// Window z centers in nm: start..end:step, start..end/count, or comma list.
    #[arg(
        long,
        value_name = "SPEC",
        default_value = "0.0..4.0:0.1",
        help_heading = "Umbrella sampling"
    )]
    pub(crate) z: String,

    /// Harmonic force constant k in kcal/mol/A^2 (GENESIS convention: U = k * (dz - z0)^2).
    #[arg(
        short = 'k',
        long,
        value_name = "K",
        default_value_t = DEFAULT_UMBRELLA_K,
        help_heading = "Umbrella sampling"
    )]
    pub(crate) k: f64,

    /// Execution mode: independent windows (default) or REUS ([REMD] type1 = RESTRAINT).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = UmbrellaModeArg::Windows,
        value_name = "MODE",
        help_heading = "Umbrella sampling"
    )]
    pub(crate) mode: UmbrellaModeArg,

    /// Dynamics ensemble for equilibration and production: nvt (default) or npt.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = UmbrellaEnsembleArg::Nvt,
        value_name = "ENSEMBLE",
        help_heading = "Umbrella sampling"
    )]
    pub(crate) ensemble: UmbrellaEnsembleArg,

    /// Replica-exchange attempt interval in MD steps (only used with --mode reus).
    #[arg(
        long,
        value_name = "STEPS",
        default_value_t = 1000,
        help_heading = "Umbrella sampling"
    )]
    pub(crate) exchange_period: usize,

    /// Restraint distance definition: distmass (COM, default) or dist (geometric center).
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = RestraintFuncArg::DistMass,
        value_name = "FUNC",
        help_heading = "Umbrella sampling"
    )]
    pub(crate) restraint_func: RestraintFuncArg,

    /// Number of blocks for MBAR error estimation.
    #[arg(
        long,
        value_name = "N",
        default_value_t = 5,
        help_heading = "Umbrella sampling"
    )]
    pub(crate) mbar_blocks: usize,

    /// PMF histogram grid spacing in Angstrom for mbar.inp.
    #[arg(
        long,
        value_name = "A",
        default_value_t = 0.5,
        help_heading = "Umbrella sampling"
    )]
    pub(crate) mbar_grid_step: f64,

    /// Output directory.
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output files"
    )]
    pub(crate) out: PathBuf,

    /// Also write per-window .pdb files.
    #[arg(long, help_heading = "Output files")]
    pub(crate) pdb: bool,

    /// Also write a shared system.psf structure file.
    #[arg(long, help_heading = "Output files")]
    pub(crate) psf: bool,

    /// Also write a VMD visualization script for window 1 (implies --psf).
    #[arg(long, help_heading = "Output files")]
    pub(crate) vmd: bool,

    #[command(flatten)]
    pub(crate) controls: ControlArgs,

    #[command(flatten)]
    pub(crate) general: GeneralArgs,
}

/// Arguments for `isolf umbrella cv`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Extract RESTR_CVS001 from GENESIS atdyn logs into cv/{}.dat for mbar_analysis",
    long_about = "\
Parse RESTR_CVS001 from completed GENESIS atdyn logs and write cv/1..N.dat for
mbar_analysis. Automatically detects independent windows (win_1..N/<log>) and
REUS trajectories (<reus-prefix>_1..N.{log,rem}, demultiplexed by parameter ID),
and supports multi-probe parent directories without requiring manifest.json."
)]
pub(crate) struct UmbrellaCvArgs {
    /// Campaign directory (or parent directory of per-probe campaigns).
    #[arg(short = 'o', long, value_name = "DIR", default_value = ".")]
    pub(crate) out: PathBuf,

    /// Equilibration steps to discard (rows with step < discard are skipped).
    #[arg(long, value_name = "STEPS", default_value_t = 0)]
    pub(crate) discard: usize,

    /// Log filename inside each win_<i>/ directory (windows mode) [default: prod.log].
    #[arg(long, value_name = "FILE", default_value = "prod.log")]
    pub(crate) log: String,

    /// Basename prefix of REUS logs and remfiles (<prefix>_<i>.{log,rem}) [default: prod_reus].
    #[arg(long, value_name = "PREFIX", default_value = "prod_reus")]
    pub(crate) reus_prefix: String,

    #[command(flatten)]
    pub(crate) general: GeneralArgs,
}

/// Arguments for top-level `isolf params`.
#[derive(clap::Args, Debug)]
#[command(
    about = "Generate or update lipid_protein.itp from parameter tables or inline overrides",
    long_about = "\
Render a GENESIS lipid_protein.itp file ([ cg_ISOLF_parameters ], auto head
[ cg_WCA_parameters ], and optional [ cg_LJ_parameters ]) in under 1 ms.

Works on any umbrella campaign directory, multi-probe parent directory,
membrane+protein directory, or standalone output file path. When --probe is
omitted and --out points to a directory, amino-acid residue types are discovered
automatically from the protein .itp files in that directory."
)]
pub(crate) struct ParamsArgs {
    /// Probe(s) or protein(s) to generate parameters for: AA code(s), 'all', or PDB/mmCIF path.
    /// If omitted, amino-acid residue types are inferred from protein .itp files in --out.
    #[arg(long, value_name = "SPEC")]
    pub(crate) probe: Option<String>,

    /// Input parameter table file ([cg_ISOLF_parameters], [cg_WCA_parameters], [cg_LJ_parameters]).
    #[arg(long, value_name = "FILE")]
    pub(crate) params: Option<PathBuf>,

    /// Override [ cg_ISOLF_parameters ] epsilon in kJ/mol.
    #[arg(long, value_name = "KJ/MOL")]
    pub(crate) eps: Option<f64>,

    /// Override [ cg_ISOLF_parameters ] sigma in nm.
    #[arg(long, value_name = "NM")]
    pub(crate) sigma: Option<f64>,

    /// Override [ cg_ISOLF_parameters ] omega in nm.
    #[arg(long, value_name = "NM")]
    pub(crate) omega: Option<f64>,

    /// Do not auto-generate Lorentz–Berthelot WCA repulsion against lipid head beads.
    #[arg(long)]
    pub(crate) no_auto_wca: bool,

    /// Target file path (e.g. lipid_protein.itp) or existing campaign / build directory.
    #[arg(short = 'o', long, value_name = "PATH", default_value = ".")]
    pub(crate) out: PathBuf,

    #[command(flatten)]
    pub(crate) general: GeneralArgs,
}

/// Arguments for `isolf lipids`.
#[derive(clap::Args, Debug)]
pub(crate) struct LipidsArgs {
    /// Disable coloured terminal output.
    #[arg(long)]
    pub(crate) no_color: bool,

    /// Use ASCII punctuation instead of Unicode glyphs.
    #[arg(long)]
    pub(crate) ascii: bool,
}
