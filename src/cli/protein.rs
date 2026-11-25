//! Standalone `isolf protein <INPUT>` coarse-graining workflow (AICG2+ / Clementi Gō).

use std::error::Error;
use std::fmt::Write as _;
use std::fs;

use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;

use isolf::force_field::ForceField;
use isolf::protein::Protein;
use isolf::protein::aicg2p::{ProteinModel, ProteinOptions, parse_charge_table, parse_selection};
use isolf::protein::params::render_aicg2p_itp;
use isolf::vmd::Layout;

use crate::report;

use super::args::{ProteinArgs, ProteinCenterArg};
use super::build::{ele_cutoff_note, write_control_files, write_vmd_scripts};
use super::parse::validate_padding;
use super::{reproduce_command, write};

pub(crate) fn run_protein(args: ProteinArgs) -> Result<(), Box<dyn Error>> {
    validate_padding(args.padding)?;
    if !args.contact_scale.is_finite() || args.contact_scale <= 0.0 {
        return Err(format!(
            "--contact-scale: '{}' must be a positive number",
            args.contact_scale
        )
        .into());
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

    let want_inp = args.want_inp();
    let want_vmd = args.want_vmd();
    let (controls, temp_note) = args.controls.prepare(want_inp, false)?;

    let prot_options = ProteinOptions {
        model: args.model.into(),
        scale: args.scale.into(),
        safe_dihedral: args.safe_dihedral.into(),
        idr_local,
        idr_nonlocal,
        charges,
        ccgo_contact_scale: args.contact_scale,
    };

    let seed = args.general.seed.unwrap_or_else(|| rand::rng().random());
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
    if let Some(note) = &temp_note {
        report::warn(note);
    }
    if want_inp && let Some(note) = ele_cutoff_note(box_size, controls.ele_cutoff) {
        report::warn(&note);
    }

    fs::create_dir_all(&args.out)?;
    let base = protein.name.clone();
    let named = |ext: &str| args.out.join(format!("{base}.{ext}"));

    write(&named("gro"), &protein.to_gro(box_size))?;
    write(&named("itp"), &protein.to_itp())?;
    write(&args.out.join("aicg2p.itp"), &render_aicg2p_itp())?;
    write(&named("top"), &protein.to_top())?;
    if args.formats.pdb {
        write(&named("pdb"), &protein.to_pdb())?;
    }
    if args.formats.psf || want_vmd {
        write(&named("psf"), &protein.to_psf())?;
    }
    if args.formats.crd {
        write(&named("crd"), &protein.to_crd())?;
    }
    if args.formats.cif {
        write(&named("cif"), &protein.to_cif())?;
    }

    if want_inp {
        write_control_files(&args.out, &base, false, box_size, &controls, &mut rng)?;
    }
    if want_vmd {
        let ff = ForceField::isolf();
        let particles = protein.to_particles(1, 1);
        write_vmd_scripts(
            &args.out,
            &base,
            Layout::Protein,
            false,
            want_inp,
            matches!(args.center, ProteinCenterArg::Origin),
            &particles,
            &ff,
            &args.vmd_opts,
            None,
            &controls,
        )?;
    }

    let log_path = named("log");
    let log_content = render_protein_log(&protein, box_size, seed, args.general.seed.is_some());
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
