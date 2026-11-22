//! Command-line builder for coarse-grained iSoLF membranes and vesicles.
//!
//! ```sh
//! isolf membrane --lipids "POPC=3,POPS=1" --count 512 --inp --vmd -o ./sim
//! isolf vesicle  --lipids POPC --radius 20 --top -o ./ves
//! isolf lipids
//! ```
//!
//! `<name>.gro` and `<name>.log` are always written. `--top` adds the topology
//! (`<name>.top` plus `isolf.itp`); `--inp` (or any `--min`/`--equil`/`--prod`
//! override) adds GENESIS control files and implies `--top`; `--vmd` (or any
//! `--vmd-*` flag) adds VMD scripts and implies `--psf`; `--gif` (or any
//! `--gif-*` flag) adds a relaxation animation; `.pdb`, `.psf`, `.crd`, and
//! `.cif` are added on request.

use std::error::Error;
use std::fs;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};

use clap::builder::styling::{AnsiColor, Styles};
use clap::{ColorChoice, CommandFactory, FromArgMatches, Parser, Subcommand};
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;

use isolf::composition::Composition;
use isolf::error::BuildError;
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

/// Build coarse-grained iSoLF lipid membranes and vesicles for GENESIS.
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

/// Options shared between `membrane` and `vesicle`.
#[derive(clap::Args, Debug)]
struct CommonBuildArgs {
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

    /// Base name for output files [default: membrane / vesicle].
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
    let (controls, temp_note) = prepare_controls(&args.common, Geometry::Membrane)?;

    let seed = args.common.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);
    let options = BuildOptions {
        padding: args.padding,
    };

    report::header();
    report::section("Build");
    let spin = report::spin("building membrane");
    let built = if args.common.want_gif() {
        Membrane::build_recorded(
            &force_field,
            Geometry::Membrane.title(),
            &upper,
            &lower,
            sizing,
            &options,
            &mut rng,
        )
    } else {
        Membrane::build(
            &force_field,
            Geometry::Membrane.title(),
            &upper,
            &lower,
            sizing,
            &options,
            &mut rng,
        )
        .map(|m| (m, Vec::new()))
    };
    spin.finish();
    let (membrane, gif_frames) = built.map_err(|e| enrich_build_error(&force_field, e))?;

    let mut notes = Vec::new();
    if let Some(n) = balance_note {
        notes.push(n);
    }
    if let Some(n) = temp_note {
        notes.push(n);
    }
    finish_build(
        membrane,
        gif_frames,
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
    let (controls, temp_note) = prepare_controls(&args.common, Geometry::Vesicle)?;

    let seed = args.common.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);
    let options = BuildOptions {
        padding: args.padding,
    };

    report::header();
    report::section("Build");
    let spin = report::spin("building vesicle");
    let built = if args.common.want_gif() {
        build_vesicle_recorded(
            &force_field,
            Geometry::Vesicle.title(),
            &outer,
            &inner,
            radius,
            &options,
            &mut rng,
        )
    } else {
        build_vesicle(
            &force_field,
            Geometry::Vesicle.title(),
            &outer,
            &inner,
            radius,
            &options,
            &mut rng,
        )
        .map(|m| (m, Vec::new()))
    };
    spin.finish();
    let (membrane, gif_frames) = built.map_err(|e| enrich_build_error(&force_field, e))?;

    let notes: Vec<String> = temp_note.into_iter().collect();
    finish_build(
        membrane,
        gif_frames,
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

#[allow(clippy::too_many_arguments)]
fn finish_build(
    mut membrane: Membrane,
    gif_frames: Vec<RelaxFrame>,
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
    report::ok(
        identity,
        &format!(
            "{} lipids {} {} beads",
            report::thousands(membrane.total_lipids()),
            report::dot(),
            report::thousands(membrane.particles.len()),
        ),
    );
    for note in notes {
        report::warn(note);
    }

    let base = common
        .name
        .clone()
        .unwrap_or_else(|| geometry.label().to_string());

    write_outputs(
        &mut membrane,
        force_field,
        common,
        geometry,
        center,
        controls,
        &base,
        &gif_frames,
        rng,
    )?;

    let out_dir = common.out.display().to_string();
    let box_size = membrane.box_size;
    let beads = membrane.particles.len();

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
}

impl Geometry {
    /// Short mode label (`"membrane"` or `"vesicle"`), also used as the default
    /// output base name when `--name` is not given.
    fn label(self) -> &'static str {
        match self {
            Geometry::Membrane => "membrane",
            Geometry::Vesicle => "vesicle",
        }
    }

    /// System title written into the structure-file headers.
    fn title(self) -> &'static str {
        match self {
            Geometry::Membrane => DEFAULT_NAME,
            Geometry::Vesicle => "CG vesicle model",
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
/// the two leaflets and counted per species (in the membrane's species order): a
/// flat membrane splits by centroid z about the mean bead z (`upper`/`lower`), a
/// vesicle by centroid radius about the mean bead radius (`outer`/`inner`).
fn leaflet_split(m: &Membrane, geometry: Geometry) -> Vec<report::Leaflet> {
    if m.particles.is_empty() {
        return Vec::new();
    }
    let vesicle = matches!(geometry, Geometry::Vesicle);
    let n = m.particles.len() as f64;
    let center = if vesicle {
        let mut c = [0.0; 3];
        for p in &m.particles {
            c[0] += p.position[0];
            c[1] += p.position[1];
            c[2] += p.position[2];
        }
        [c[0] / n, c[1] / n, c[2] / n]
    } else {
        [0.0; 3]
    };
    let scalar = |pos: [f64; 3]| {
        if vesicle {
            let dx = pos[0] - center[0];
            let dy = pos[1] - center[1];
            let dz = pos[2] - center[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        } else {
            pos[2]
        }
    };

    let threshold = m.particles.iter().map(|p| scalar(p.position)).sum::<f64>() / n;
    let n_species = m.lipid_counts.len();
    let mut far = vec![0usize; n_species];
    let mut near = vec![0usize; n_species];
    let mut species_idx = 0usize;
    for mol in m.particles.chunk_by(|a, b| a.molecule_id == b.molecule_id) {
        let res = mol[0].residue.as_str();
        if species_idx >= n_species || m.lipid_counts[species_idx].0 != res {
            if let Some(idx) = m.lipid_counts.iter().position(|(s, _)| s == res) {
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
        m.lipid_counts
            .iter()
            .zip(counts)
            .filter(|&(_, &c)| c > 0)
            .map(|((species, _), &c)| (species.clone(), c))
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
    membrane: &mut Membrane,
    force_field: &ForceField,
    common: &CommonBuildArgs,
    geometry: Geometry,
    center: CenterArg,
    controls: &ControlSpecs,
    base: &str,
    gif_frames: &[RelaxFrame],
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&common.out)?;
    let named = |extension: &str| common.out.join(format!("{base}.{extension}"));

    membrane.center_in_box();
    let box_size = membrane.box_size;
    if matches!(center, CenterArg::Origin) {
        membrane.translate([-box_size[0] / 2.0, -box_size[1] / 2.0, -box_size[2] / 2.0]);
    }
    let membrane: &Membrane = membrane;

    let want_inp = common.want_inp();
    let want_vmd = common.want_vmd();
    let want_gif = common.want_gif();

    let formats: [StructureFormat; 5] = [
        (true, "gro", Membrane::to_gro),
        (common.pdb, "pdb", Membrane::to_pdb),
        (common.psf || want_vmd, "psf", Membrane::to_psf),
        (common.crd, "crd", Membrane::to_crd),
        (common.cif, "cif", Membrane::to_cif),
    ];
    for (enabled, ext, render) in formats {
        if enabled {
            write(&named(ext), &render(membrane))?;
        }
    }

    if common.top || want_inp {
        write(&named("top"), &membrane.to_top())?;
        write(&common.out.join("isolf.itp"), &force_field.to_itp())?;
    }

    if want_gif && !gif_frames.is_empty() {
        let bounds = match geometry {
            Geometry::Membrane => [0.0, 0.0, box_size[0], box_size[1]],
            Geometry::Vesicle => [-2.0, -2.0, 2.0, 2.0],
        };
        let opts = anim::GifOptions {
            mode: match common.gif_mode.unwrap_or(GifModeArg::Point) {
                GifModeArg::Point => anim::GifMode::Point,
                GifModeArg::Density => anim::GifMode::Density,
            },
            scale: common.gif_scale.unwrap_or(1.0),
            fps: common.gif_fps.unwrap_or(14),
        };
        let path = named("gif");
        anim::write_gif(&path, gif_frames, bounds, &opts)?;
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
            membrane,
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
    membrane: &Membrane,
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
            &membrane.particles,
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

        // Discovery subcommand.
        let cli = Cli::try_parse_from(["isolf", "lipids"]).unwrap();
        assert!(matches!(cli.command, Command::Lipids(_)));
    }
}
