use clap::Parser;
use colored::Colorize;
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
}

fn main() {
    let args = Args::parse();
    let output_path = PathBuf::from(&args.output);

    if output_path.exists() {
        eprintln!(
            "{} Output path '{}' already exists",
            "Error:".red().bold(),
            output_path.to_string_lossy().blue().bold()
        );
        std::process::exit(1);
    }

    if let Err(e) = fs::create_dir_all(&output_path) {
        eprintln!("{} {}", "Error:".red().bold(), e);
        std::process::exit(1);
    }

    let equ_input_file = isolf::InputFileBuilder::default()
        .input_grotop("./membrane.top")
        .input_grocrd("./membrane.gro")
        .output_dcd(isolf::Output::new("./equilibration.dcd", 100))
        .output_rst(isolf::Output::new("./equilibration.rst", 10000))
        .num_steps(10000)
        .output_ene_period(100)
        .update_nb_period(10)
        .remove_tr_period(10)
        .ensemble(isolf::Ensemble::npt(303.15, 0.0, 0.01, 0.01))
        .boundary(isolf::Boundary::pbc_with_box_size(255.821, 255.821, 200.0))
        .build()
        .unwrap_or_else(|e| {
            eprintln!("{} {}", "Error:".red().bold(), e);
            std::process::exit(1);
        });
    let equ_input_path = output_path.join("equilibration.inp");
    if let Err(e) = fs::write(&equ_input_path, equ_input_file.to_string()) {
        eprintln!("{} {}", "Error:".red().bold(), e);
        std::process::exit(1);
    }

    let pro_input_file = isolf::InputFileBuilder::default()
        .input_grotop("./membrane.top")
        .input_grocrd("./membrane.gro")
        .input_rst("./equilibration.rst")
        .output_dcd(isolf::Output::new("./production.dcd", 100))
        .output_rst(isolf::Output::new("./production.rst", 10000))
        .num_steps(10000)
        .output_ene_period(100)
        .update_nb_period(10)
        .remove_tr_period(10)
        .ensemble(isolf::Ensemble::npt(303.15, 0.0, 0.01, 0.01))
        .boundary(isolf::Boundary::pbc())
        .build()
        .unwrap_or_else(|e| {
            eprintln!("{} {}", "Error:".red().bold(), e);
            std::process::exit(1);
        });
    let pro_input_path = output_path.join("production.inp");
    if let Err(e) = fs::write(&pro_input_path, pro_input_file.to_string()) {
        eprintln!("{} {}", "Error:".red().bold(), e);
        std::process::exit(1);
    }
}
