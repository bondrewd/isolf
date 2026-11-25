//! Command-line dispatch, terminal styling, and shared CLI helpers.

mod args;
mod build;
mod lipids;
mod parse;
mod protein;
mod umbrella;

use std::error::Error;
use std::fs;
use std::io::IsTerminal;
use std::path::Path;

use clap::builder::styling::{AnsiColor, Styles};
use clap::{ColorChoice, CommandFactory, FromArgMatches};

use isolf::force_field::ForceField;

use crate::{report, update};

use args::{Cli, Command, UmbrellaCommand};

/// Upper bound for a randomly drawn Langevin thermostat seed.
pub(crate) const MAX_ISEED: u64 = 999_999;

/// Colour scheme for the `--help` output.
fn help_styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Green.on_default().bold())
        .usage(AnsiColor::Green.on_default().bold())
        .literal(AnsiColor::Cyan.on_default().bold())
        .placeholder(AnsiColor::Cyan.on_default())
}

pub(crate) fn main() {
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
            a.common.general.ascii,
            a.common.general.quiet,
            a.common.general.verbose,
            a.common.general.no_color,
        ),
        Command::Vesicle(a) => (
            a.common.general.ascii,
            a.common.general.quiet,
            a.common.general.verbose,
            a.common.general.no_color,
        ),
        Command::Tubule(a) => (
            a.common.general.ascii,
            a.common.general.quiet,
            a.common.general.verbose,
            a.common.general.no_color,
        ),
        Command::Protein(a) => (
            a.general.ascii,
            a.general.quiet,
            a.general.verbose,
            a.general.no_color,
        ),
        Command::Umbrella(UmbrellaCommand::Setup(a)) => (
            a.general.ascii,
            a.general.quiet,
            a.general.verbose,
            a.general.no_color,
        ),
        Command::Umbrella(UmbrellaCommand::Cv(a)) => (
            a.general.ascii,
            a.general.quiet,
            a.general.verbose,
            a.general.no_color,
        ),
        Command::Params(a) => (
            a.general.ascii,
            a.general.quiet,
            a.general.verbose,
            a.general.no_color,
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
        Command::Membrane(args) => build::run_membrane(args),
        Command::Vesicle(args) => build::run_vesicle(args),
        Command::Tubule(args) => build::run_tubule(args),
        Command::Protein(args) => protein::run_protein(args),
        Command::Umbrella(UmbrellaCommand::Setup(args)) => umbrella::run_umbrella_setup(*args),
        Command::Umbrella(UmbrellaCommand::Cv(args)) => umbrella::run_umbrella_cv(args),
        Command::Params(args) => umbrella::run_params(args),
        Command::Lipids(_) => {
            lipids::print_lipids(&ForceField::isolf(), color, unicode);
            Ok(())
        }
        Command::Update => update::run(),
        Command::Uninstall => update::uninstall(),
    }
}

/// Reconstruct the invocation as a re-runnable, shell-quoted line, with the seed
/// pinned in (appended when the user did not pass `--seed`).
pub(crate) fn reproduce_command(seed: u64, seed_given: bool) -> String {
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

pub(crate) fn write(path: &Path, contents: &str) -> std::io::Result<()> {
    fs::write(path, contents)?;
    report::record(path);
    Ok(())
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::args::{PhaseArg, ProteinModelArg, UmbrellaModeArg};
    use super::*;

    #[test]
    fn clap_definition_passes_debug_assert() {
        Cli::command().debug_assert();
    }

    #[test]
    fn cli_subcommands_and_auto_enabled_flags() {
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
            "--lipid-protein",
            "lp.txt",
            "--ele-cutoff",
            "35",
        ])
        .unwrap();
        let Command::Membrane(args) = cli.command else {
            panic!("expected membrane subcommand");
        };
        assert!(args.common.want_inp());
        assert!(args.common.want_vmd());
        assert!(args.common.want_gif());
        assert!(args.common.lipid_protein.is_some());
        assert_eq!(args.common.controls.ele_cutoff, Some(35.0));

        // --lipids conflicts with --upper.
        assert!(
            Cli::try_parse_from([
                "isolf", "membrane", "-L", "POPC", "--upper", "DOPC", "-n", "64"
            ])
            .is_err()
        );

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

        // Umbrella setup + cv + top-level params subcommands.
        let cli = Cli::try_parse_from([
            "isolf",
            "umbrella",
            "setup",
            "--probe",
            "walp23.pdb,tilt=90",
            "-L",
            "POPC",
            "-n",
            "64",
            "--z",
            "0.0..4.0:0.1",
            "-k",
            "0.5",
            "--mode",
            "reus",
        ])
        .unwrap();
        let Command::Umbrella(UmbrellaCommand::Setup(umb)) = cli.command else {
            panic!("expected umbrella setup subcommand");
        };
        assert_eq!(umb.probe, "walp23.pdb,tilt=90");
        assert_eq!(umb.mode, UmbrellaModeArg::Reus);
        assert_eq!(umb.k, 0.5);

        // --bilayer conflicts with -L / -n.
        assert!(
            Cli::try_parse_from([
                "isolf",
                "umbrella",
                "setup",
                "--probe",
                "TRP",
                "--bilayer",
                "eq.gro",
                "-L",
                "POPC",
            ])
            .is_err()
        );

        let cli = Cli::try_parse_from([
            "isolf",
            "umbrella",
            "cv",
            "--discard",
            "1000",
            "-o",
            "us_trp",
        ])
        .unwrap();
        let Command::Umbrella(UmbrellaCommand::Cv(cv)) = cli.command else {
            panic!("expected umbrella cv subcommand");
        };
        assert_eq!(cv.discard, 1000);

        let cli = Cli::try_parse_from([
            "isolf",
            "params",
            "--params",
            "trial.txt",
            "--probe",
            "TRP",
            "-o",
            "us_trp/lipid_protein.itp",
        ])
        .unwrap();
        assert!(matches!(cli.command, Command::Params(_)));

        let cli = Cli::try_parse_from(["isolf", "lipids"]).unwrap();
        assert!(matches!(cli.command, Command::Lipids(_)));
    }
}
