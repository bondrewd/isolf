//! Command-line membrane builder.
//!
//! Builds a simulation-ready iSoLF coarse-grained bilayer from options passed
//! directly on the command line and writes the output files:
//!
//! ```sh
//! isolf --upper POPC=1,DOPS=1 --lipids-per-leaflet 512 --out ./sim --pdb --psf --inp
//! ```
//!
//! `<name>.gro` is always written. `--top` adds the topology (`<name>.top` plus
//! `isolf.itp`); `--inp` adds the GENESIS control files and implies `--top`;
//! `.pdb`, `.psf`, `.crd`, `.cif`, and `.vmd` are added on request.

use std::error::Error;
use std::fs;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};

use clap::builder::styling::{AnsiColor, Styles};
use clap::{ArgGroup, ColorChoice, CommandFactory, FromArgMatches, Parser, Subcommand};
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;

use isolf::composition::Composition;
use isolf::force_field::ForceField;
use isolf::membrane::{
    BuildOptions, DEFAULT_NAME, DEFAULT_PADDING, DEFAULT_TEMPERATURE, Membrane, RelaxFrame, Sizing,
};
use isolf::vesicle::{VesicleRadius, build_vesicle, build_vesicle_recorded};
use isolf::vmd::{Layout, VmdOptions, VmdSource};

mod anim;
mod report;
mod runlog;
mod update;

/// Where the finished system is centred (`--center`). The periodic box is the
/// same either way; only the coordinates move.
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq)]
enum CenterArg {
    /// Centred on the coordinate origin `(0, 0, 0)` (default). The VMD scripts
    /// translate the system into the box for display, so the view is unchanged.
    #[default]
    Origin,
    /// Centred in the box `[0, box]`, matching the drawn periodic box.
    Box,
}

/// Scope for `--vmd-color-by`: how lipid beads are coloured. Maps to
/// [`isolf::vmd::LipidColoring`].
#[derive(clap::ValueEnum, Clone, Copy, Debug)]
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
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq)]
enum GifModeArg {
    /// One disc per lipid; clear for small to medium systems.
    #[default]
    Point,
    /// Local-density heatmap; better for large, dense systems.
    Density,
}

/// Maintenance subcommands. With none given, isolf builds a structure from the
/// options below.
#[derive(Subcommand, Debug)]
enum Command {
    /// Update isolf in place to the latest release, then exit.
    Update(Maintenance),
    /// Uninstall isolf, removing the binary from the system, then exit.
    Uninstall(Maintenance),
}

/// Flags shared by the maintenance subcommands: just a manual `--help`. The parent
/// disables its auto help flag (to reposition it under "Extra options") and clap
/// propagates that down, so without this `isolf update --help` would not work.
#[derive(clap::Args, Debug)]
struct Maintenance {
    /// Print help
    #[arg(short = 'h', long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// Build coarse-grained lipid membranes and vesicles for GENESIS.
#[derive(Parser, Debug)]
// Help/version are added explicitly (at the bottom of the struct) so the "Extra
// options" group renders last in --help instead of clap's default of first.
#[command(
    name = "isolf",
    version,
    about,
    disable_help_flag = true,
    disable_version_flag = true,
    disable_help_subcommand = true
)]
#[command(long_about = "\
Build coarse-grained starting structures for GENESIS molecular dynamics.

The flags you give select one of two build modes:
  membrane   --upper <lipids> --lipids-per-leaflet <N>  (or --membrane)
  vesicle    --upper <lipids> --vesicle <R>

The structure (.gro) goes to --out. Add --top for the topology (.top + .itp),
--inp for GENESIS control files (implies --top), and --vmd for VMD scripts.")]
// A geometry and a composition are validated together in run(); both groups are
// optional so the cross-checks (a geometry needs a composition, and vice versa)
// can report a helpful message instead of clap's terse one.
#[command(group = ArgGroup::new("geometry").required(false).multiple(true)
    .args(["lipids_per_leaflet", "membrane", "vesicle"]))]
#[command(group = ArgGroup::new("composition").required(false).multiple(true)
    .args(["upper", "lower"]))]
struct Args {
    // ---- Lipid composition ----
    /// Upper leaflet (the outer shell for a vesicle) as NAME=WEIGHT pairs.
    ///
    /// Examples:
    /// --upper "POPC"                    # Upper leaflet contains POPC
    /// --upper "POPC=1,DOPS=1"           # Upper leaflet contains 50% of POPC and 50% of DOPS
    /// --upper "POPC=1/2,DOPS=1/2"       # Using fractions is also possible
    /// --upper "POPC=50%,DOPS=50%"       # Using percentages is also possible
    /// --upper "POPC=1" --upper "DOPS=1" # Adding one by one is also possible
    #[arg(
        short = 'u',
        long,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Lipid composition"
    )]
    upper: Vec<String>,

    /// Lower leaflet (the inner shell for a vesicle) as NAME=WEIGHT pairs.
    ///
    /// Examples:
    /// --lower "POPC"                    # Lower leaflet contains POPC
    /// --lower "POPC=1,DOPS=1"           # Lower leaflet contains 50% of POPC and 50% of DOPS
    /// --lower "POPC=1/2,DOPS=1/2"       # Using fractions is also possible
    /// --lower "POPC=50%,DOPS=50%"       # Using percentages is also possible
    /// --lower "POPC=1" --lower "DOPS=1" # Adding one by one is also possible
    #[arg(
        short = 'l',
        long,
        value_name = "SPEC",
        verbatim_doc_comment,
        help_heading = "Lipid composition"
    )]
    lower: Vec<String>,

    // ---- Geometry ----
    /// Number of lipids per leaflet in membranes.
    ///
    /// Examples:
    /// --lipids-per-leaflet 1000                                    # Add 1000 lipids in each leaflet
    /// --lipids-per-leaflet "up=1200,lo=800"                        # Add 1200 lipids to the upper leaflet and 800 to the lower
    /// --lipids-per-leaflet "up=1200" --lipids-per-leaflet "lo=800" # Adding one by one is also possible
    #[arg(
        long,
        value_name = "SPEC",
        conflicts_with = "vesicle",
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    lipids_per_leaflet: Vec<String>,

    /// Flat-membrane size in nm as "x=10,y=20".
    ///
    /// Examples:
    /// --membrane 10                       # Make a membrane of 10 nm square (x = y = 10)
    /// --membrane "x=10"                   # Make a membrane of x = 10 nm and make x = y
    /// --membrane "x=10,y=20"              # Make a membrane of x = 10 nm by y = 20 nm
    /// --membrane "x=10" --membrane "y=20" # Adding one by one is also possible
    #[arg(
        long,
        value_name = "SPEC",
        conflicts_with = "vesicle",
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    membrane: Vec<String>,

    /// Vesicle radius in nm as "ro=20" (outer) or "ri=15" (inner).
    ///
    /// Examples:
    /// --vesicle 20                        # Make a vesicle with an outer radius of 20 nm and determine the inner radius automatically
    /// --vesicle "ro=20"                   # Make a vesicle with an outer radius of 20 nm and determine the inner radius automatically
    /// --vesicle "ri=15"                   # Make a vesicle with an inner radius of 15 nm and determine the outer radius automatically
    /// --vesicle "ri=15,ro=20"             # Make a vesicle with an inner radius of 15 nm and an outer radius of 20 nm
    /// --vesicle "ri=15" --vesicle "ro=20" # Adding one by one is also possible
    #[arg(
        long,
        value_name = "SPEC",
        conflicts_with_all = ["lipids_per_leaflet", "membrane"],
        verbatim_doc_comment,
        help_heading = "Geometry"
    )]
    vesicle: Vec<String>,

    /// System center position. Possible values are "origin" (0,0,0) and "box" (Lx/2,Ly/2,Lz/2)
    #[arg(long, value_enum, ignore_case = true, default_value_t = CenterArg::Origin, hide_possible_values = true, value_name = "WHERE", help_heading = "Geometry")]
    center: CenterArg,

    /// Padding in nm along z for membranes and on all sides for vesicles
    #[arg(long, value_name = "P", default_value_t = DEFAULT_PADDING, help_heading = "Geometry")]
    padding: f64,

    /// When both leaflet counts are set and their areas differ, the lighter leaflet is grown to match.
    ///
    /// When setting this flag, the program doesn't grow the smaller leaflet where there is a size difference.
    #[arg(long, verbatim_doc_comment, help_heading = "Geometry")]
    no_balance: bool,

    // ---- Output files ----
    /// Output directory for saving the output
    #[arg(
        short = 'o',
        long,
        value_name = "DIR",
        default_value = ".",
        help_heading = "Output files"
    )]
    out: PathBuf,

    /// Base name for the structure/topology files.
    ///
    /// [default: membrane/vesicle depending on the geometry]
    #[arg(long, value_name = "NAME", help_heading = "Output files")]
    name: Option<String>,

    /// Also write a topology ".top" file and the force-field ".itp" file.
    ///
    /// When using the "--inp" flag, this flag is set automatically as well since
    /// the topology is necessary for running the MD simulation.
    #[arg(long, help_heading = "Output files")]
    top: bool,

    /// Also write a ".pdb" coordinate file
    #[arg(long, help_heading = "Output files")]
    pdb: bool,

    /// Also write a CHARMM ".psf" structure file
    #[arg(long, help_heading = "Output files")]
    psf: bool,

    /// Also write a CHARMM ".crd" coordinate file
    #[arg(long, help_heading = "Output files")]
    crd: bool,

    /// Also write an mmCIF ".cif" coordinate file
    #[arg(long, help_heading = "Output files")]
    cif: bool,

    /// Also write an animation of the leaflet relaxation.
    #[arg(long, help_heading = "Output files")]
    gif: bool,

    /// Also write GENESIS ".inp" control files.
    ///
    /// Three files are generated depending on the geometry.
    ///
    /// For membranes:
    /// - min.inp: Energy minimization
    /// - npt.inp: Equilibration in the NPT ensemble
    /// - pro.inp: Production run in the NPT ensemble
    ///
    /// For vesicles:
    /// - min.inp: Energy minimization
    /// - nvt.inp: Equilibration in the NVT ensemble
    /// - pro.inp: Production run in the NVT ensemble
    #[arg(long, verbatim_doc_comment, help_heading = "Output files")]
    inp: bool,

    /// Also write VMD ".vmd" scripts for visualization.
    ///
    /// To execute the scripts, use the "-e" flag in VMD.
    ///
    /// Example:
    /// vmd -e membrane.vmd
    #[arg(long, verbatim_doc_comment, help_heading = "Output files")]
    vmd: bool,

    // ---- General options ----
    /// System temperature in K for the GENESIS ".inp" files.
    ///
    /// The temperature of each simulation can be overwritten by specifying it
    /// with any of the following flags:
    /// --nvt-temperature
    /// --npt-temperature
    /// --pro-temperature
    #[arg(short = 't', long, value_name = "K", default_value_t = DEFAULT_TEMPERATURE, verbatim_doc_comment, help_heading = "General options")]
    temperature: f64,

    /// Seed value for the RNG
    ///
    /// A random number is used if not assigned
    #[arg(long, value_name = "SEED", help_heading = "General options")]
    seed: Option<u64>,

    #[command(flatten)]
    min: MinControl,

    #[command(flatten)]
    nvt: NvtControl,

    #[command(flatten)]
    npt: NptControl,

    #[command(flatten)]
    pro: ProControl,

    #[command(flatten)]
    viz: VmdControl,

    // ---- Relaxation visualization (--gif) ----
    /// Visualization mode.
    ///
    /// How `--gif` renders: a "point" scatter or a "density" heatmap.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        default_value_t = GifModeArg::Point,
        value_name = "MODE",
        help_heading = "Relaxation visualization (--gif)"
    )]
    gif_mode: GifModeArg,

    /// Resolution scale (e.g. 2 doubles it).
    ///
    /// Scale the `--gif` resolution (e.g. 2 doubles it).
    #[arg(
        long,
        value_name = "FACTOR",
        default_value_t = 1.0,
        help_heading = "Relaxation visualization (--gif)"
    )]
    gif_scale: f64,

    /// Playback speed in frames per second.
    ///
    /// Playback speed of the `--gif`, in frames per second.
    #[arg(
        long,
        value_name = "FPS",
        default_value_t = 14,
        help_heading = "Relaxation visualization (--gif)"
    )]
    gif_fps: u16,

    // The standard options, declared last (and the auto help/version disabled
    // above) so the "Extra options" group renders at the bottom of --help.
    /// Print only the final summary and any warnings.
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        help_heading = "Extra options"
    )]
    quiet: bool,

    /// Print each file written and per-phase detail.
    #[arg(short = 'v', long, help_heading = "Extra options")]
    verbose: bool,

    /// Disable coloured terminal output.
    #[arg(long, help_heading = "Extra options")]
    no_color: bool,

    /// Use ASCII status markers instead of Unicode glyphs.
    #[arg(long, help_heading = "Extra options")]
    ascii: bool,

    /// Maintenance commands: `isolf update` and `isolf uninstall`.
    #[command(subcommand)]
    command: Option<Command>,

    /// Print help (see more with '--help').
    #[arg(short = 'h', long, action = clap::ArgAction::Help, help_heading = "Extra options")]
    help: Option<bool>,

    /// Print version.
    #[arg(short = 'V', long, action = clap::ArgAction::Version, help_heading = "Extra options")]
    version: Option<bool>,
}

/// Minimization control-file parameters (used with `--inp`).
#[derive(clap::Args, Debug)]
#[command(next_help_heading = "Minimization control file (--inp)")]
struct MinControl {
    /// Number of minimization steps.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    min_num_steps: usize,
    /// Energy output period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    min_eneout_period: usize,
    /// Coordinate output period.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    min_crdout_period: usize,
    /// Restart output period.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    min_rstout_period: usize,
    /// Neighbour-list update period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    min_nbupdate_period: usize,
}

/// NVT control-file parameters (used with `--inp`).
#[derive(clap::Args, Debug)]
#[command(next_help_heading = "NVT control file (--inp)")]
struct NvtControl {
    /// Number of NVT steps.
    #[arg(long, value_name = "N", default_value_t = 20000)]
    nvt_num_steps: usize,
    /// NVT time step.
    #[arg(long, value_name = "DT", default_value_t = 0.01)]
    nvt_time_step: f64,
    /// NVT temperature in K (default to --temperature).
    #[arg(long, value_name = "K")]
    nvt_temperature: Option<f64>,
    /// Energy output period.
    #[arg(long, value_name = "N", default_value_t = 100)]
    nvt_eneout_period: usize,
    /// Coordinate output period.
    #[arg(long, value_name = "N", default_value_t = 100)]
    nvt_crdout_period: usize,
    /// Restart output period.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    nvt_rstout_period: usize,
    /// Translational-motion removal period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    nvt_stoptr_period: usize,
    /// Neighbour-list update period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    nvt_nbupdate_period: usize,
    /// Langevin thermostat seed (random if unset).
    #[arg(long, value_name = "SEED")]
    nvt_seed: Option<u64>,
}

/// NPT control-file parameters (used with `--inp`).
#[derive(clap::Args, Debug)]
#[command(next_help_heading = "NPT control file (--inp)")]
struct NptControl {
    /// Number of NPT steps.
    #[arg(long, value_name = "N", default_value_t = 20000)]
    npt_num_steps: usize,
    /// NPT time step.
    #[arg(long, value_name = "DT", default_value_t = 0.01)]
    npt_time_step: f64,
    /// NPT temperature in K (default to --temperature).
    #[arg(long, value_name = "K")]
    npt_temperature: Option<f64>,
    /// Energy output period.
    #[arg(long, value_name = "N", default_value_t = 100)]
    npt_eneout_period: usize,
    /// Coordinate output period.
    #[arg(long, value_name = "N", default_value_t = 100)]
    npt_crdout_period: usize,
    /// Restart output period.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    npt_rstout_period: usize,
    /// Translational-motion removal period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    npt_stoptr_period: usize,
    /// Neighbour-list update period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    npt_nbupdate_period: usize,
    /// Langevin thermostat seed (random if unset).
    #[arg(long, value_name = "SEED")]
    npt_seed: Option<u64>,
}

/// Production control-file parameters (used with `--inp`): the run continues from
/// the equilibration restart under the same ensemble (NPT for a membrane, NVT for
/// a vesicle).
#[derive(clap::Args, Debug)]
#[command(next_help_heading = "Production control file (--inp)")]
struct ProControl {
    /// Number of production steps.
    #[arg(long, value_name = "N", default_value_t = 500000)]
    pro_num_steps: usize,
    /// Production time step.
    #[arg(long, value_name = "DT", default_value_t = 0.02)]
    pro_time_step: f64,
    /// Production temperature in K (default to --temperature).
    #[arg(long, value_name = "K")]
    pro_temperature: Option<f64>,
    /// Energy output period.
    #[arg(long, value_name = "N", default_value_t = 500)]
    pro_eneout_period: usize,
    /// Coordinate output period.
    #[arg(long, value_name = "N", default_value_t = 500)]
    pro_crdout_period: usize,
    /// Restart output period.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    pro_rstout_period: usize,
    /// Translational-motion removal period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    pro_stoptr_period: usize,
    /// Neighbour-list update period.
    #[arg(long, value_name = "N", default_value_t = 10)]
    pro_nbupdate_period: usize,
    /// Langevin thermostat seed (random if unset).
    #[arg(long, value_name = "SEED")]
    pro_seed: Option<u64>,
}

/// VMD visualization options (used with `--vmd`).
#[derive(clap::Args, Debug)]
#[command(next_help_heading = "VMD visualization (--vmd)")]
struct VmdControl {
    /// Use a smooth QuickSurf surface instead of VDW spheres.
    #[arg(long)]
    vmd_surface: bool,
    /// Slice the system with a clipping plane to reveal the cross-section.
    #[arg(long)]
    vmd_cutaway: bool,
    /// Append a Tachyon ray-tracing render command to each script.
    #[arg(long)]
    vmd_render: bool,
    /// Background colour.
    ///
    /// A VMD colour name.
    #[arg(long, value_name = "COLOR", default_value = "white")]
    vmd_background: String,
    /// Representation material.
    ///
    /// A VMD material name, e.g. AOChalky or Goodsell.
    #[arg(long, value_name = "NAME", default_value = "AOChalky")]
    vmd_material: String,
    /// Trajectory smoothing window, in frames.
    #[arg(long, value_name = "N", default_value_t = 5)]
    vmd_smooth: usize,
    /// Colouring scheme.
    ///
    /// `head` tints the head bead by species (rest by role), `bead` colours every
    /// bead by role, `whole` colours the whole lipid by species.
    #[arg(
        long,
        value_enum,
        ignore_case = true,
        value_name = "SCOPE",
        default_value = "head"
    )]
    vmd_color_by: LipidColorScope,
}

impl NvtControl {
    /// Build the GENESIS dynamics parameters from these NVT flags, drawing a
    /// random thermostat seed when one was not given explicitly.
    fn dynamics(&self, default_temperature: f64, rng: &mut StdRng) -> isolf::inp::Dynamics {
        isolf::inp::Dynamics {
            temperature: self.nvt_temperature.unwrap_or(default_temperature),
            num_steps: self.nvt_num_steps,
            time_step: self.nvt_time_step,
            eneout_period: self.nvt_eneout_period,
            crdout_period: self.nvt_crdout_period,
            rstout_period: self.nvt_rstout_period,
            stoptr_period: self.nvt_stoptr_period,
            nbupdate_period: self.nvt_nbupdate_period,
            seed: self
                .nvt_seed
                .unwrap_or_else(|| rng.random_range(1..=MAX_ISEED)),
        }
    }
}

impl NptControl {
    /// Build the GENESIS dynamics parameters from these NPT flags.
    fn dynamics(&self, default_temperature: f64, rng: &mut StdRng) -> isolf::inp::Dynamics {
        isolf::inp::Dynamics {
            temperature: self.npt_temperature.unwrap_or(default_temperature),
            num_steps: self.npt_num_steps,
            time_step: self.npt_time_step,
            eneout_period: self.npt_eneout_period,
            crdout_period: self.npt_crdout_period,
            rstout_period: self.npt_rstout_period,
            stoptr_period: self.npt_stoptr_period,
            nbupdate_period: self.npt_nbupdate_period,
            seed: self
                .npt_seed
                .unwrap_or_else(|| rng.random_range(1..=MAX_ISEED)),
        }
    }
}

impl ProControl {
    /// Build the GENESIS dynamics parameters from these production flags.
    fn dynamics(&self, default_temperature: f64, rng: &mut StdRng) -> isolf::inp::Dynamics {
        isolf::inp::Dynamics {
            temperature: self.pro_temperature.unwrap_or(default_temperature),
            num_steps: self.pro_num_steps,
            time_step: self.pro_time_step,
            eneout_period: self.pro_eneout_period,
            crdout_period: self.pro_crdout_period,
            rstout_period: self.pro_rstout_period,
            stoptr_period: self.pro_stoptr_period,
            nbupdate_period: self.pro_nbupdate_period,
            seed: self
                .pro_seed
                .unwrap_or_else(|| rng.random_range(1..=MAX_ISEED)),
        }
    }
}

/// Upper bound for a randomly drawn Langevin thermostat seed.
const MAX_ISEED: u64 = 999_999;

/// One optional structure-file format: whether it was requested, its extension,
/// and how to render it from the membrane.
type StructureFormat = (bool, &'static str, fn(&Membrane) -> String);

/// Colour scheme for the `--help` output.
fn help_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Green.on_default().bold())
        .usage(AnsiColor::Green.on_default().bold())
        .literal(AnsiColor::Cyan.on_default().bold())
        .placeholder(AnsiColor::Cyan.on_default())
}

fn main() {
    // Detected before parsing so it can also colour `--help` and error output.
    let no_color =
        std::env::args().any(|arg| arg == "--no-color") || std::env::var_os("NO_COLOR").is_some();
    let color = if no_color {
        ColorChoice::Never
    } else {
        ColorChoice::Auto
    };
    let command = Args::command().color(color).styles(help_styles());
    let args = Args::from_arg_matches(&command.get_matches()).unwrap_or_else(|e| e.exit());
    let color_on = !no_color && std::io::stdout().is_terminal();
    report::init(color_on, !args.ascii, args.quiet, args.verbose);

    if let Err(error) = run(args) {
        report::error(&error.to_string());
        std::process::exit(1);
    }
}

fn run(args: Args) -> Result<(), Box<dyn Error>> {
    if let Some(command) = &args.command {
        return match command {
            Command::Update(_) => update::run(),
            Command::Uninstall(_) => update::uninstall(),
        };
    }
    report::header();
    let force_field = ForceField::isolf();
    // Always materialise a seed (random when none is given) so the run log can
    // record a command that reproduces this build exactly.
    let seed = args.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);

    let membrane_counts = parse_lipids_per_leaflet(&args.lipids_per_leaflet)?;
    let membrane_box = parse_membrane(&args.membrane)?;
    let vesicle = parse_vesicle(&args.vesicle)?;
    let want_membrane = membrane_counts.is_some() || membrane_box.is_some() || vesicle.is_some();

    // Combine the (repeatable, comma-separated) leaflet specs into one string each.
    let upper_input = join_specs(&args.upper);
    let lower_input = join_specs(&args.lower);
    let upper_given = !upper_input.is_empty();
    let lower_given = !lower_input.is_empty();

    // A composition with no size can't build a membrane; refuse rather than
    // silently dropping the lipids.
    if (upper_given || lower_given) && !want_membrane {
        return Err("a lipid composition (--upper/--lower) needs a size:\n  \
             membrane:  --lipids-per-leaflet <N>  (or --membrane <L>)\n  \
             vesicle:   --vesicle <R>"
            .into());
    }
    if !want_membrane {
        return Err("nothing to build. Choose a mode:\n  \
             membrane:  --upper <lipids> --lipids-per-leaflet <N>  (or --membrane <L>)\n  \
             vesicle:   --upper <lipids> --vesicle <R>"
            .into());
    }

    // The build mode follows from the flags: a vesicle if a radius was given,
    // otherwise a flat membrane.
    let geometry = if vesicle.is_some() {
        Geometry::Vesicle
    } else {
        Geometry::Membrane
    };
    report::section("Build");

    // A leaflet may be omitted; it then copies the other.
    let upper_opt = upper_given.then_some(&upper_input);
    let lower_opt = lower_given.then_some(&lower_input);
    let upper_spec = upper_opt
        .or(lower_opt)
        .ok_or("a membrane needs a composition (--upper or --lower)")?;
    let lower_spec = lower_opt.or(upper_opt).unwrap_or(upper_spec);
    let upper = Composition::parse(upper_spec)?;
    let lower = Composition::parse(lower_spec)?;
    let options = BuildOptions {
        padding: args.padding,
    };
    let kind = if geometry == Geometry::Vesicle {
        "vesicle"
    } else {
        "membrane"
    };
    let spin = report::spin(&format!("building {kind}"));
    let (membrane, note, gif_frames) = match geometry {
        Geometry::Vesicle => {
            let radius = vesicle.expect("a vesicle geometry has a radius");
            let (membrane, frames) = if args.gif {
                build_vesicle_recorded(
                    &force_field,
                    geometry.title(),
                    &upper,
                    &lower,
                    radius,
                    &options,
                    &mut rng,
                )?
            } else {
                let membrane = build_vesicle(
                    &force_field,
                    geometry.title(),
                    &upper,
                    &lower,
                    radius,
                    &options,
                    &mut rng,
                )?;
                (membrane, Vec::new())
            };
            (membrane, None, frames)
        }
        Geometry::Membrane => {
            let (sizing, note) = resolve_membrane_sizing(
                &args,
                &force_field,
                &upper,
                &lower,
                membrane_counts,
                membrane_box,
            )?;
            let (membrane, frames) = if args.gif {
                Membrane::build_recorded(
                    &force_field,
                    geometry.title(),
                    &upper,
                    &lower,
                    sizing,
                    &options,
                    &mut rng,
                )?
            } else {
                let membrane = Membrane::build(
                    &force_field,
                    geometry.title(),
                    &upper,
                    &lower,
                    sizing,
                    &options,
                    &mut rng,
                )?;
                (membrane, Vec::new())
            };
            (membrane, note, frames)
        }
    };
    spin.finish();
    report::ok(
        kind,
        &format!(
            "{} lipids {} {} beads",
            report::thousands(membrane.total_lipids()),
            report::dot(),
            report::thousands(membrane.particles.len()),
        ),
    );
    if let Some(note) = note {
        report::warn(&note);
    }

    // Output file base name: `--name`, else the build-mode default. Shared by the
    // written files and the run log below.
    let base = args
        .name
        .clone()
        .unwrap_or_else(|| geometry.default_base_name().to_string());

    write_outputs(
        &membrane,
        &force_field,
        &args,
        geometry,
        &base,
        &gif_frames,
        &mut rng,
    )?;

    // Box, particle count, and mode label, shared by the run log and the summary.
    let out_dir = args.out.display().to_string();
    let box_size = membrane.box_size;
    let beads = membrane.particles.len();
    let identity = mode_label(geometry);

    // The always-on run log, written next to the output and listed within it.
    let files: Vec<String> = report::recorded_files()
        .iter()
        .filter_map(|p| p.file_name().map(|f| f.to_string_lossy().into_owned()))
        .collect();
    let command = reproduce_command(seed, args.seed.is_some());
    let log_path = args.out.join(format!("{base}.log"));
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

    // The grouped Write section (now including the log) and the closing summary.
    report::write_section(&out_dir);
    let next = args.inp.then(|| {
        let equil = if geometry.is_npt() {
            "npt.inp"
        } else {
            "nvt.inp"
        };
        // atdyn resolves the .inp's relative include paths against the working
        // directory, so it must run inside the output dir. Prefix a `cd` only when
        // a real `--out` was given (skip the pointless `cd .` for the default).
        let cd = if out_dir == "." {
            String::new()
        } else {
            format!("cd {out_dir} && ")
        };
        format!("run in GENESIS: {cd}atdyn min.inp (then {equil}, pro.inp)")
    });
    let leaflets = leaflet_split(&membrane, geometry);
    report::summary(&report::Summary {
        identity,
        box_nm: box_size,
        lipids: membrane.lipid_counts.as_slice(),
        leaflets: Some(&leaflets),
        beads,
        next,
    });
    Ok(())
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
}

impl Geometry {
    /// System title written into the structure-file headers.
    fn title(self) -> &'static str {
        match self {
            Geometry::Membrane => DEFAULT_NAME,
            Geometry::Vesicle => "CG vesicle model",
        }
    }

    /// Default base name for the output files when `--name` is not given.
    fn default_base_name(self) -> &'static str {
        match self {
            Geometry::Membrane => "membrane",
            Geometry::Vesicle => "vesicle",
        }
    }

    /// VMD camera layout (a vesicle gets the perspective view).
    fn layout(self) -> Layout {
        match self {
            Geometry::Membrane => Layout::Membrane,
            Geometry::Vesicle => Layout::Vesicle,
        }
    }

    /// Whether the production/equilibration ensemble is NPT (a flat membrane) as
    /// opposed to NVT (a vesicle, at fixed volume).
    fn is_npt(self) -> bool {
        self == Geometry::Membrane
    }
}

/// Parse a single positive numeric dimension (nm), trimming whitespace.
fn parse_dim(value: &str) -> Result<f64, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v > 0.0)
        .ok_or_else(|| format!("'{value}' is not a positive number of nm").into())
}

/// Parse a positive lipid count, trimming whitespace.
fn parse_count(value: &str) -> Result<usize, Box<dyn Error>> {
    let value = value.trim();
    value
        .parse::<usize>()
        .ok()
        .filter(|n| *n > 0)
        .ok_or_else(|| format!("'{value}' is not a positive lipid count").into())
}

/// Parse the repeatable `--lipids-per-leaflet` occurrences into per-leaflet lipid
/// counts `(upper, lower)`, or `None` if the flag was absent. Each occurrence is a
/// comma-separated list of `up=<N>`, `lo=<N>`, or a bare number (sets `up`); an
/// unspecified leaflet copies the other, so `1000` and `up=1000` both mean 1000
/// lipids in each leaflet.
fn parse_lipids_per_leaflet(specs: &[String]) -> Result<Option<(usize, usize)>, Box<dyn Error>> {
    let (mut up, mut lo) = (None, None);
    let mut seen = false;
    for token in specs.iter().flat_map(|s| s.split(',')) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        seen = true;
        match token.split_once('=') {
            Some((key, value)) => match key.trim().to_ascii_lowercase().as_str() {
                "up" | "upper" => up = Some(parse_count(value)?),
                "lo" | "lower" => lo = Some(parse_count(value)?),
                other => {
                    return Err(format!(
                        "--lipids-per-leaflet: unknown key '{other}' (use up or lo)"
                    )
                    .into());
                }
            },
            None => up = Some(parse_count(token)?),
        }
    }
    if !seen {
        return Ok(None);
    }
    let counts = match (up, lo) {
        (Some(up), Some(lo)) => (up, lo),
        (Some(n), None) | (None, Some(n)) => (n, n),
        (None, None) => {
            return Err(
                "--lipids-per-leaflet needs a value, e.g. `1000` or `up=1000,lo=2000`".into(),
            );
        }
    };
    Ok(Some(counts))
}

/// Parse the repeatable `--membrane` occurrences into a flat-membrane box
/// `(x, y)` in nm, or `None` if the flag was absent. Each occurrence is a
/// comma-separated list of `x=<nm>`, `y=<nm>`, or a bare number (sets `x`); an
/// unspecified side copies the other, so `10` and `x=10` both mean a 10 nm square.
fn parse_membrane(specs: &[String]) -> Result<Option<(f64, f64)>, Box<dyn Error>> {
    let (mut x, mut y) = (None, None);
    let mut seen = false;
    for token in specs.iter().flat_map(|s| s.split(',')) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        seen = true;
        match token.split_once('=') {
            Some((key, value)) => match key.trim().to_ascii_lowercase().as_str() {
                "x" => x = Some(parse_dim(value)?),
                "y" => y = Some(parse_dim(value)?),
                other => {
                    return Err(format!("--membrane: unknown key '{other}' (use x or y)").into());
                }
            },
            None => x = Some(parse_dim(token)?),
        }
    }
    if !seen {
        return Ok(None);
    }
    let dims = match (x, y) {
        (Some(x), Some(y)) => (x, y),
        (Some(s), None) | (None, Some(s)) => (s, s),
        (None, None) => {
            return Err("--membrane needs a value, e.g. `10` or `x=10,y=20`".into());
        }
    };
    Ok(Some(dims))
}

/// Parse the repeatable `--vesicle` occurrences into a [`VesicleRadius`], or
/// `None` if the flag was absent. Each occurrence is a comma-separated list of
/// `ri=<nm>` (inner/lumen radius), `ro=<nm>` (outer radius), or a bare number (sets
/// the outer radius). With one radius the other is derived from the bilayer
/// thickness; with both, [`build_vesicle`] checks the gap holds the bilayer.
fn parse_vesicle(specs: &[String]) -> Result<Option<VesicleRadius>, Box<dyn Error>> {
    let (mut inner, mut outer) = (None, None);
    let mut seen = false;
    for token in specs.iter().flat_map(|s| s.split(',')) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        seen = true;
        match token.split_once('=') {
            Some((key, value)) => match key.trim().to_ascii_lowercase().as_str() {
                "ri" => inner = Some(parse_dim(value)?),
                "ro" => outer = Some(parse_dim(value)?),
                other => {
                    return Err(format!("--vesicle: unknown key '{other}' (use ro or ri)").into());
                }
            },
            None => outer = Some(parse_dim(token)?),
        }
    }
    if !seen {
        return Ok(None);
    }
    let radius = match (inner, outer) {
        (Some(inner), Some(outer)) => VesicleRadius::Both { inner, outer },
        (Some(inner), None) => VesicleRadius::Inner(inner),
        (None, Some(outer)) => VesicleRadius::Outer(outer),
        (None, None) => {
            return Err("--vesicle needs a value, e.g. `20` or `ri=15,ro=20`".into());
        }
    };
    Ok(Some(radius))
}

/// Combine the repeatable, comma-separated occurrences of a leaflet flag
/// (`--upper POPC=1 --upper DPPC=1` or `--upper "POPC=1,DPPC=1"`) into a single
/// `NAME=WEIGHT,…` string, dropping blank occurrences. The combined string is then
/// parsed by [`Composition::parse`].
fn join_specs(specs: &[String]) -> String {
    specs
        .iter()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join(",")
}

/// Map the count/box options to a [`Sizing`]:
/// counts only → minimum square box; box only → fill it; both → pack into the box.
fn resolve_sizing(
    counts: Option<(usize, usize)>,
    box_xy: Option<(f64, f64)>,
) -> Result<Sizing, Box<dyn Error>> {
    match (counts, box_xy) {
        (Some((upper, lower)), None) => Ok(Sizing::Count { upper, lower }),
        (None, Some((x, y))) => Ok(Sizing::Box { x, y }),
        (Some((upper, lower)), Some((x, y))) => Ok(Sizing::CountInBox { upper, lower, x, y }),
        // A flat-membrane geometry always carries a count or a box.
        (None, None) => unreachable!("a membrane geometry has a count or a box"),
    }
}

/// Resolve the flat-membrane sizing, matching the two leaflets' areas. Returns
/// the sizing and an optional note for the user. Balancing applies only to count
/// mode; an explicit box already fills each leaflet at its own target density.
fn resolve_membrane_sizing(
    args: &Args,
    force_field: &ForceField,
    upper: &Composition,
    lower: &Composition,
    counts: Option<(usize, usize)>,
    box_xy: Option<(f64, f64)>,
) -> Result<(Sizing, Option<String>), Box<dyn Error>> {
    let (Some((upper_count, lower_count)), None) = (counts, box_xy) else {
        return Ok((resolve_sizing(counts, box_xy)?, None));
    };
    let upper_apl = isolf::membrane::area_per_lipid(force_field, upper)?;
    let lower_apl = isolf::membrane::area_per_lipid(force_field, lower)?;
    let (upper, lower, note) = balanced_counts(
        upper_apl,
        lower_apl,
        upper_count,
        lower_count,
        !args.no_balance,
    );
    Ok((Sizing::Count { upper, lower }, note))
}

/// Resolve the per-leaflet counts so the leaflets start at matched area (the same
/// lateral density), from each leaflet's area per lipid (nm²).
///
/// If the two leaflets' areas differ by more than 2%, the counts are kept and a
/// note returned (the lighter leaflet would start under-dense, a possible pore),
/// unless `balance` grows the lighter leaflet to match. Returns (upper, lower, note).
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
        // Grow the lighter leaflet to the heavier leaflet's area.
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
             (a possible pore). Set {key}={to_match} in --lipids-per-leaflet (from {count_now}), \
             or drop --no-balance to grow it to match automatically"
        );
        (upper, lower, Some(note))
    }
}

/// A short label for the build mode, shown in the closing summary.
fn mode_label(geometry: Geometry) -> &'static str {
    match geometry {
        Geometry::Membrane => "membrane",
        Geometry::Vesicle => "vesicle",
    }
}

/// Per-leaflet lipid composition. Each lipid molecule is classified into one of
/// the two leaflets and counted per species (in the membrane's species order): a
/// flat membrane splits by centroid z about the mean bead z (`upper`/`lower`), a
/// vesicle by centroid radius about the mean bead radius (`outer`/`inner`).
fn leaflet_split(m: &Membrane, geometry: Geometry) -> Vec<report::Leaflet> {
    use std::collections::HashMap;
    if m.particles.is_empty() {
        return Vec::new();
    }
    let vesicle = matches!(geometry, Geometry::Vesicle);
    let n = m.particles.len() as f64;
    // The scalar each bead is classified by: distance from the lipid centre (a
    // vesicle's radius) or its z (a flat membrane).
    let center = {
        let mut c = [0.0; 3];
        for p in &m.particles {
            for (acc, v) in c.iter_mut().zip(p.position) {
                *acc += v / n;
            }
        }
        c
    };
    let scalar = |pos: [f64; 3]| {
        if vesicle {
            let d = [pos[0] - center[0], pos[1] - center[1], pos[2] - center[2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        } else {
            pos[2]
        }
    };
    // Midplane threshold and each molecule's species + mean scalar.
    let threshold = m.particles.iter().map(|p| scalar(p.position)).sum::<f64>() / n;
    let mut molecules: HashMap<usize, (String, f64, usize)> = HashMap::new();
    for p in &m.particles {
        let entry = molecules
            .entry(p.molecule_id)
            .or_insert_with(|| (p.residue.clone(), 0.0, 0));
        entry.1 += scalar(p.position);
        entry.2 += 1;
    }
    let mut far: HashMap<&str, usize> = HashMap::new();
    let mut near: HashMap<&str, usize> = HashMap::new();
    for (species, sum, count) in molecules.values() {
        let side = if *sum / *count as f64 >= threshold {
            &mut far
        } else {
            &mut near
        };
        *side.entry(species.as_str()).or_insert(0) += 1;
    }
    let in_species_order = |counts: &HashMap<&str, usize>| -> Vec<(String, usize)> {
        m.lipid_counts
            .iter()
            .filter_map(|(species, _)| {
                let c = counts.get(species.as_str()).copied().unwrap_or(0);
                (c > 0).then(|| (species.clone(), c))
            })
            .collect()
    };
    let (far_label, near_label) = if vesicle {
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
    membrane: &Membrane,
    force_field: &ForceField,
    args: &Args,
    geometry: Geometry,
    base: &str,
    gif_frames: &[RelaxFrame],
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&args.out)?;
    let named = |extension: &str| args.out.join(format!("{base}.{extension}"));

    // The flat bilayer is built centred on the origin, which draws in a corner of
    // the box, so centre it in the box for the viewer (a no-op for a vesicle,
    // already centred).
    let centered = membrane.centered_in_box();

    // The periodic box is fixed by the build; `--center` only moves the
    // coordinates within it.
    let box_size = centered.box_size;

    // `--center origin` shifts the box-centred build to the coordinate origin
    // (minus half the box). The box — and so the `.gro` box line — is unchanged;
    // the VMD scripts translate the system back into the box, so the view matches
    // the box-centred default. Box mode leaves the coordinates exactly as built.
    let shifted = matches!(args.center, CenterArg::Origin).then(|| {
        let mut m = centered.clone();
        m.translate([-box_size[0] / 2.0, -box_size[1] / 2.0, -box_size[2] / 2.0]);
        m
    });
    let membrane = shifted.as_ref().unwrap_or(&centered);

    // The `.gro` is always written.
    write(&named("gro"), &membrane.to_gro())?;

    // Topology files (`.top` plus the force-field `.itp`) with --top, or implied
    // by --inp, since a run needs them.
    if args.top || args.inp {
        write(&named("top"), &membrane.to_top())?;
        write(&args.out.join("isolf.itp"), &force_field.to_itp())?;
    }

    // Optional structure formats. --vmd also needs the .psf.
    let formats: [StructureFormat; 4] = [
        (args.pdb, "pdb", Membrane::to_pdb),
        (args.psf || args.vmd, "psf", Membrane::to_psf),
        (args.crd, "crd", Membrane::to_crd),
        (args.cif, "cif", Membrane::to_cif),
    ];
    for (enabled, ext, to_membrane) in formats {
        if enabled {
            write(&named(ext), &to_membrane(membrane))?;
        }
    }

    // The relaxation animation for --gif. A membrane records box (x, y) positions;
    // a vesicle records a Lambert azimuthal equal-area map in the disc of radius 2.
    if args.gif && !gif_frames.is_empty() {
        let bounds = match geometry {
            Geometry::Membrane => [0.0, 0.0, box_size[0], box_size[1]],
            Geometry::Vesicle => [-2.0, -2.0, 2.0, 2.0],
        };
        let opts = anim::GifOptions {
            mode: match args.gif_mode {
                GifModeArg::Point => anim::GifMode::Point,
                GifModeArg::Density => anim::GifMode::Density,
            },
            scale: args.gif_scale,
            fps: args.gif_fps,
        };
        let path = named("gif");
        anim::write_gif(&path, gif_frames, bounds, &opts)?;
        report::record(&path);
    }

    if args.inp {
        write_control_files(&args.out, base, geometry, box_size, args, rng)?;
    }
    if args.vmd {
        write_vmd_scripts(&args.out, base, geometry, membrane, force_field, args)?;
    }
    Ok(())
}

/// Write the GENESIS control files: minimize, then equilibrate and run
/// production under the build mode's ensemble — NPT for a flat membrane (which
/// relaxes the area), NVT for a vesicle (fixed volume). Production continues from
/// the equilibration restart.
fn write_control_files(
    out: &Path,
    base: &str,
    geometry: Geometry,
    box_size: [f64; 3],
    args: &Args,
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    let temperature = args.temperature;
    let minimization = isolf::inp::Minimization {
        box_size,
        temperature,
        num_steps: args.min.min_num_steps,
        eneout_period: args.min.min_eneout_period,
        crdout_period: args.min.min_crdout_period,
        rstout_period: args.min.min_rstout_period,
        nbupdate_period: args.min.min_nbupdate_period,
    };
    write(
        &out.join("min.inp"),
        &isolf::inp::minimization(base, &minimization)?,
    )?;

    let production = args.pro.dynamics(temperature, rng);
    if geometry.is_npt() {
        let npt = args.npt.dynamics(temperature, rng);
        write(
            &out.join("npt.inp"),
            &isolf::inp::npt(base, "min.rst", "npt", &npt)?,
        )?;
        write(
            &out.join("pro.inp"),
            &isolf::inp::npt(base, "npt.rst", "pro", &production)?,
        )?;
    } else {
        let nvt = args.nvt.dynamics(temperature, rng);
        write(
            &out.join("nvt.inp"),
            &isolf::inp::nvt(base, "min.rst", "nvt", &nvt)?,
        )?;
        write(
            &out.join("pro.inp"),
            &isolf::inp::nvt(base, "nvt.rst", "pro", &production)?,
        )?;
    }
    Ok(())
}

/// Write one VMD visualization script per simulation stage, each loading that
/// stage's trajectory. The frame count is the stage's nsteps / coordinate-output
/// period.
fn write_vmd_scripts(
    out: &Path,
    base: &str,
    geometry: Geometry,
    membrane: &Membrane,
    force_field: &ForceField,
    args: &Args,
) -> Result<(), Box<dyn Error>> {
    let options = VmdOptions {
        surface: args.viz.vmd_surface,
        cutaway: args.viz.vmd_cutaway,
        render: args.viz.vmd_render,
        background: &args.viz.vmd_background,
        material: &args.viz.vmd_material,
        smooth: args.viz.vmd_smooth,
        lipid_coloring: match args.viz.vmd_color_by {
            LipidColorScope::Head => isolf::vmd::LipidColoring::Head,
            LipidColorScope::Bead => isolf::vmd::LipidColoring::Role,
            LipidColorScope::Whole => isolf::vmd::LipidColoring::Whole,
        },
        recenter_to_box: matches!(args.center, CenterArg::Origin),
    };
    let layout = geometry.layout();
    let script = |source| {
        isolf::vmd::vmd_script(
            force_field,
            &membrane.particles,
            base,
            layout,
            source,
            &options,
        )
    };

    // Always: a script for the initial built structure (`<base>.vmd`).
    write(
        &out.join(format!("{base}.vmd")),
        &script(VmdSource::Structure),
    )?;

    // Only with --inp: one script per simulation stage's trajectory. The
    // equilibration stage is NPT for a flat membrane, NVT otherwise.
    if args.inp {
        let (equil, equil_steps, equil_crdout) = if geometry.is_npt() {
            ("npt", args.npt.npt_num_steps, args.npt.npt_crdout_period)
        } else {
            ("nvt", args.nvt.nvt_num_steps, args.nvt.nvt_crdout_period)
        };
        let stages = [
            ("min", args.min.min_num_steps, args.min.min_crdout_period),
            (equil, equil_steps, equil_crdout),
            ("pro", args.pro.pro_num_steps, args.pro.pro_crdout_period),
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
        // Symmetric (equal area per lipid, equal counts): unchanged, no note.
        let (u, l, note) = balanced_counts(0.64, 0.64, 200, 200, false);
        assert_eq!((u, l), (200, 200));
        assert!(note.is_none());

        // Mismatch, no balance: keep counts, suggest the fix via --lipids-per-leaflet.
        let (u, l, note) = balanced_counts(0.64, 0.64, 200, 150, false);
        assert_eq!((u, l), (200, 150));
        assert!(note.unwrap().contains("lo=200 in --lipids-per-leaflet"));

        // Mismatch, balance on: grow the lighter leaflet to match.
        let (u, l, note) = balanced_counts(0.64, 0.64, 200, 150, true);
        assert_eq!((u, l), (200, 200));
        assert!(note.unwrap().contains("balanced"));
    }
}
