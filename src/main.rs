use anyhow::{Result, anyhow};
use clap::Parser;
use colored::Colorize;
use isolf::itp::{force_field::ForceField, itp_file::ItpFile};
use std::fs;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "isolf")]
#[command(about = "Tool for initializing iSoLF CG MD simulations")]
#[command(version)]
struct Args {
    /// Output directory path
    #[arg(short)]
    #[arg(long)]
    #[arg(value_name = "PATH")]
    #[arg(default_value_t = String::from("sim"))]
    output: String,

    /// Simulation temperature
    #[arg(long)]
    #[arg(value_name = "TEMPERATURE")]
    #[arg(default_value_t = 303.15)]
    temperature: f64,

    /// Equilibration run steps
    #[arg(long)]
    #[arg(value_name = "STEPS")]
    #[arg(default_value_t = 10_000)]
    equilibration_steps: u64,

    /// Equilibration rst output period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 1000)]
    equilibration_rst_period: u64,

    /// Equilibration dcd output period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 1000)]
    equilibration_dcd_period: u64,

    /// Equilibration ene output period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 100)]
    equilibration_ene_period: u64,

    /// Equilibration neighbor list update period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 10)]
    equilibration_nb_period: u64,

    /// Equilibration translation and rotation removal period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 10)]
    equilibration_tr_period: u64,

    /// Production run steps
    #[arg(long)]
    #[arg(value_name = "STEPS")]
    #[arg(default_value_t = 100_000)]
    production_steps: u64,

    /// Production rst output period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 1000)]
    production_rst_period: u64,

    /// Production dcd output period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 1000)]
    production_dcd_period: u64,

    /// Production ene output period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 100)]
    production_ene_period: u64,

    /// Production neighbor list update period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 10)]
    production_nb_period: u64,

    /// Production translation and rotation removal period
    #[arg(long)]
    #[arg(value_name = "PERIOD")]
    #[arg(default_value_t = 10)]
    production_tr_period: u64,
}

fn run() -> Result<()> {
    let args = Args::parse();
    let output_path = PathBuf::from(&args.output);

    if output_path.exists() {
        return Err(anyhow!(
            "Output path '{}' already exists",
            output_path.to_string_lossy().blue().bold()
        ));
    }

    fs::create_dir_all(&output_path)?;

    let equ_input_file = isolf::inp::builder::InputFileBuilder::default()
        .input_grotop("./membrane.top")
        .input_grocrd("./membrane.gro")
        .output_dcd("./equilibration.dcd", args.equilibration_dcd_period)
        .output_rst("./equilibration.rst", args.equilibration_rst_period)
        .solvent_temperature(args.temperature)
        .num_steps(args.equilibration_steps)
        .output_ene_period(args.equilibration_ene_period)
        .update_nb_period(args.equilibration_nb_period)
        .remove_tr_period(args.equilibration_tr_period)
        .npt(args.temperature, 0.0, 0.01, 0.01)
        .pbc_with_box_size(255.821, 255.821, 200.0)
        .build()?;
    let equ_input_path = output_path.join("equilibration.inp");
    fs::write(&equ_input_path, equ_input_file.to_string())?;

    let pro_input_file = isolf::inp::builder::InputFileBuilder::default()
        .input_grotop("./membrane.top")
        .input_grocrd("./membrane.gro")
        .input_rst("./equilibration.rst")
        .output_dcd("./production.dcd", args.production_dcd_period)
        .output_rst("./production.rst", args.equilibration_rst_period)
        .solvent_temperature(args.temperature)
        .num_steps(args.production_steps)
        .output_ene_period(args.production_ene_period)
        .update_nb_period(args.production_nb_period)
        .remove_tr_period(args.production_tr_period)
        .npt(args.temperature, 0.0, 0.01, 0.01)
        .pbc()
        .build()?;
    let pro_input_path = output_path.join("production.inp");
    fs::write(&pro_input_path, pro_input_file.to_string())?;

    let isolf_ff: ForceField = serde_json::from_str(include_str!("../data/ff.json"))?;
    let isolf_itp = ItpFile::try_from(isolf_ff)?;
    println!("{}", isolf_itp);

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{} {}", "Error:".red().bold(), e);
        std::process::exit(1);
    }
}
