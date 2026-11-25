//! CLI runners for `isolf umbrella setup`, `isolf umbrella cv`, and `isolf params`.

use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rayon::prelude::*;

use isolf::embed::EmbeddedSystem;
use isolf::force_field::{ForceField, LipidProteinParams, ProbeParamOverrides};
use isolf::inp::{MbarConfig, Minimization, ReusConfig, ZRestraint};
use isolf::membrane::{BuildOptions, DEFAULT_NAME, Membrane};
use isolf::protein::params::{canonical_aa, render_aicg2p_itp};
use isolf::umbrella::{self, CvSource, UmbrellaManifest, UmbrellaProbe};
use isolf::vmd::{Layout, VmdOptions, VmdSource};

use crate::report;

use super::args::{
    ParamsArgs, UmbrellaCvArgs, UmbrellaEnsembleArg, UmbrellaModeArg, UmbrellaSetupArgs,
};
use super::build::{ele_cutoff_note, resolve_membrane_sizing};
use super::lipids::enrich_build_error;
use super::parse::{
    ControlSpecs, parse_box_spec, parse_count_spec, parse_umbrella_probes, resolve_compositions,
};
use super::{MAX_ISEED, write};

/// Validated inputs for `isolf umbrella setup`, resolved up front before any
/// files are written.
struct ValidatedSetup {
    probes: Vec<UmbrellaProbe>,
    resolved_params: Vec<LipidProteinParams>,
    z_windows: Vec<f64>,
    controls: ControlSpecs,
    min_box_z: f64,
}

fn validate_setup(
    args: &UmbrellaSetupArgs,
    force_field: &ForceField,
) -> Result<ValidatedSetup, Box<dyn Error>> {
    let probes = parse_umbrella_probes(&args.probe, args.charge)?;
    if !args.k.is_finite() || args.k <= 0.0 {
        return Err(format!("--k: '{}' must be a positive force constant", args.k).into());
    }
    if !args.box_z.is_finite() || args.box_z <= 0.0 {
        return Err(format!("--box-z: '{}' must be a positive number of nm", args.box_z).into());
    }
    if args.mbar_blocks == 0 {
        return Err("--mbar-blocks must be at least 1".into());
    }
    if !args.mbar_grid_step.is_finite() || args.mbar_grid_step <= 0.0 {
        return Err(format!(
            "--mbar-grid-step: '{}' must be positive",
            args.mbar_grid_step
        )
        .into());
    }
    let z_windows = umbrella::parse_z_windows(&args.z)?;
    if args.mode == UmbrellaModeArg::Reus && z_windows.len() < 2 {
        return Err("--mode reus requires at least 2 umbrella windows".into());
    }
    if args.mode == UmbrellaModeArg::Reus && args.exchange_period == 0 {
        return Err("--exchange-period must be positive for --mode reus".into());
    }

    let is_npt = matches!(args.ensemble, UmbrellaEnsembleArg::Npt);
    let (controls, _) = args.controls.prepare(true, is_npt)?;
    if args.mode == UmbrellaModeArg::Reus {
        if !controls
            .equil
            .num_steps
            .is_multiple_of(args.exchange_period)
        {
            return Err(format!(
                "--equil steps ({}) must be a multiple of --exchange-period ({})",
                controls.equil.num_steps, args.exchange_period
            )
            .into());
        }
        if !controls.prod.num_steps.is_multiple_of(args.exchange_period) {
            return Err(format!(
                "--prod steps ({}) must be a multiple of --exchange-period ({})",
                controls.prod.num_steps, args.exchange_period
            )
            .into());
        }
    }

    let base_params = match &args.params {
        Some(path) => {
            let content = fs::read_to_string(path)
                .map_err(|e| format!("--lipid-protein: failed to read {}: {e}", path.display()))?;
            Some(LipidProteinParams::parse(&content)?)
        }
        None => None,
    };
    let overrides = ProbeParamOverrides {
        epsilon: args.eps,
        sigma: args.sigma,
        omega: args.omega,
        auto_head_wca: !args.no_auto_wca,
    };

    let mut resolved_params = Vec::with_capacity(probes.len());
    for probe in &probes {
        resolved_params.push(LipidProteinParams::for_proteins(
            force_field,
            &[&probe.protein],
            base_params.as_ref(),
            overrides,
        )?);
    }

    let max_abs_z = z_windows.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let max_probe_half_z = probes
        .iter()
        .map(|p| {
            let (lo, hi) = p.protein.bounding_box();
            0.5 * (hi[2] - lo[2])
        })
        .fold(0.0_f64, f64::max);
    let min_box_z = args.box_z.max(2.0 * (max_abs_z + max_probe_half_z + 5.0));

    Ok(ValidatedSetup {
        probes,
        resolved_params,
        z_windows,
        controls,
        min_box_z,
    })
}

fn load_or_build_bilayer(
    args: &UmbrellaSetupArgs,
    force_field: &ForceField,
    min_box_z: f64,
    rng: &mut StdRng,
) -> Result<Membrane, Box<dyn Error>> {
    let mut bilayer = if let Some(gro_path) = &args.bilayer_gro {
        let content = fs::read_to_string(gro_path)
            .map_err(|e| format!("--bilayer: failed to read {}: {e}", gro_path.display()))?;
        let parsed = isolf::gro::parse_bilayer(&content, force_field)?;
        report::header();
        report::section("Build");
        parsed
    } else {
        let (upper, lower) = resolve_compositions(
            &args.bilayer.lipids,
            &args.bilayer.upper,
            &args.bilayer.lower,
            "--upper",
            "--lower",
        )?;
        let counts = parse_count_spec(&args.bilayer.count)?;
        let box_xy = parse_box_spec(&args.bilayer.box_size)?;
        let counts = if counts.is_none() && box_xy.is_none() {
            Some((64, 64))
        } else {
            counts
        };
        let (sizing, _) =
            resolve_membrane_sizing(false, force_field, &upper, &lower, counts, box_xy)
                .map_err(|e| enrich_build_error(force_field, e))?;
        let build_opts = BuildOptions {
            padding: 0.5 * (min_box_z - 4.0).max(4.0),
            ..BuildOptions::default()
        };
        report::header();
        report::section("Build");
        let spin = report::spin("building reference bilayer");
        let built = Membrane::build(
            force_field,
            DEFAULT_NAME,
            &upper,
            &lower,
            sizing,
            &build_opts,
            rng,
        )
        .map_err(|e| enrich_build_error(force_field, e))?;
        spin.finish();
        built
    };

    umbrella::center_bilayer_at_origin(&mut bilayer, min_box_z);
    Ok(bilayer)
}

struct RenderedWindow {
    id: usize,
    z_a: f64,
    gro: String,
    pdb: Option<String>,
}

/// Build and format all windows in parallel while preserving exact deterministic
/// order by window index. Returns window 1's [`EmbeddedSystem`] for shared
/// `.top`/`.psf`/`.vmd` rendering plus all pre-rendered window coordinate files.
fn build_and_render_windows(
    bilayer: &Membrane,
    probe: &UmbrellaProbe,
    z_windows: &[f64],
    label: &str,
    want_pdb: bool,
    force_field: &ForceField,
) -> (EmbeddedSystem, Vec<RenderedWindow>) {
    let probe_name = probe.protein.name.as_str();
    let built: Vec<(EmbeddedSystem, RenderedWindow)> = z_windows
        .par_iter()
        .enumerate()
        .map(|(idx, &z_nm)| {
            let id = idx + 1;
            let title = format!("{probe_name} {label} {id} z={z_nm:.3}nm");
            let sys =
                umbrella::build_window_system(title, bilayer, &probe.protein, z_nm, force_field);
            let gro = sys.to_gro();
            let pdb = want_pdb.then(|| sys.to_pdb());
            let rw = RenderedWindow {
                id,
                z_a: z_nm * 10.0,
                gro,
                pdb,
            };
            (sys, rw)
        })
        .collect();

    let mut iter = built.into_iter();
    let (first_sys, first_rw) = iter.next().expect("z_windows is validated to be non-empty");
    let mut rendered = Vec::with_capacity(z_windows.len());
    rendered.push(first_rw);
    rendered.extend(iter.map(|(_, rw)| rw));
    (first_sys, rendered)
}

struct CampaignContext<'a> {
    args: &'a UmbrellaSetupArgs,
    force_field: &'a ForceField,
    bilayer: &'a Membrane,
    controls: &'a ControlSpecs,
    min_params: &'a Minimization,
    z_windows: &'a [f64],
    z_angstrom: &'a [f64],
    k_vec: &'a [f64],
    group2: &'a str,
    isolf_itp: &'a str,
    aicg2p_itp: &'a str,
    grid_min_a: f64,
    grid_max_a: f64,
    grid_points: usize,
}

fn write_probe_campaign(
    ctx: &CampaignContext<'_>,
    probe: &UmbrellaProbe,
    lp_params: &LipidProteinParams,
    campaign_dir: &Path,
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    let probe_name = probe.protein.name.as_str();
    let group1 = probe.selection();
    let func_str = ctx.args.restraint_func.as_genesis_str();
    let ensemble = ctx.args.ensemble.to_ensemble();
    let is_reus = ctx.args.mode == UmbrellaModeArg::Reus;
    let mode_label = if is_reus { "reus" } else { "windows" };

    fs::create_dir_all(campaign_dir)?;
    write(&campaign_dir.join("aicg2p.itp"), ctx.aicg2p_itp)?;
    write(&campaign_dir.join("isolf.itp"), ctx.isolf_itp)?;
    write(
        &campaign_dir.join("lipid_protein.itp"),
        &lp_params.render_itp(),
    )?;
    write(
        &campaign_dir.join(format!("{probe_name}.itp")),
        &probe.protein.to_itp(),
    )?;

    let (first_sys, windows) = build_and_render_windows(
        ctx.bilayer,
        probe,
        ctx.z_windows,
        if is_reus { "replica" } else { "window" },
        ctx.args.pdb,
        ctx.force_field,
    );

    write(
        &campaign_dir.join("system.top"),
        &first_sys.to_top_with_lipid_protein(true),
    )?;
    if ctx.args.psf || ctx.args.vmd {
        write(&campaign_dir.join("system.psf"), &first_sys.to_psf())?;
    }
    if ctx.args.vmd {
        let vmd_gro = if is_reus {
            "./win_1.gro"
        } else {
            "./win_1/system.gro"
        };
        let vmd_opts = VmdOptions {
            surface: false,
            cutaway: false,
            render: false,
            background: "white",
            material: "AOChalky",
            smooth: 5,
            lipid_coloring: isolf::vmd::LipidColoring::Head,
            recenter_to_box: true,
        };
        let vmd = isolf::vmd::vmd_script(
            ctx.force_field,
            &first_sys.particles,
            "system",
            Layout::Membrane,
            VmdSource::StructureAt { gro: vmd_gro },
            &vmd_opts,
        );
        write(&campaign_dir.join("system.vmd"), &vmd)?;
    }

    if is_reus {
        for win in &windows {
            let gro_name = format!("win_{}.gro", win.id);
            write(&campaign_dir.join(&gro_name), &win.gro)?;
            if let Some(pdb) = &win.pdb {
                write(&campaign_dir.join(format!("win_{}.pdb", win.id)), pdb)?;
            }
            let restraint = ZRestraint {
                group1: &group1,
                group2: ctx.group2,
                function: func_str,
                direction: "Z2",
                constant: ctx.args.k,
                reference_angstrom: win.z_a,
            };
            write(
                &campaign_dir.join(format!("min_{}.inp", win.id)),
                &isolf::inp::minimization_with_restraint(
                    "./system.top",
                    &format!("./{gro_name}"),
                    &format!("min_{}", win.id),
                    ctx.min_params,
                    Some(&restraint),
                )?,
            )?;
        }
        let eq_seed = rng.random_range(1..=MAX_ISEED);
        let eq_remd_seed = rng.random_range(1..=MAX_ISEED);
        let pr_seed = rng.random_range(1..=MAX_ISEED);
        let pr_remd_seed = rng.random_range(1..=MAX_ISEED);
        let eq_reus = ReusConfig {
            exchange_period: ctx.args.exchange_period,
            iseed: eq_remd_seed,
            group1: &group1,
            group2: ctx.group2,
            function: func_str,
            direction: "Z2",
            constants: ctx.k_vec,
            references_angstrom: ctx.z_angstrom,
        };
        let pr_reus = ReusConfig {
            iseed: pr_remd_seed,
            ..eq_reus
        };
        write(
            &campaign_dir.join("equil_reus.inp"),
            &isolf::inp::reus(
                ensemble,
                "./system.top",
                "./win_{}.gro",
                "./min_{}.rst",
                "equil_reus",
                &ctx.controls.equil.to_params(
                    ctx.controls.temperature,
                    eq_seed,
                    ctx.bilayer.box_size,
                    ctx.controls.ele_cutoff,
                ),
                &eq_reus,
            )?,
        )?;
        write(
            &campaign_dir.join("prod_reus.inp"),
            &isolf::inp::reus(
                ensemble,
                "./system.top",
                "./win_{}.gro",
                "./equil_reus_{}.rst",
                "prod_reus",
                &ctx.controls.prod.to_params(
                    ctx.controls.temperature,
                    pr_seed,
                    ctx.bilayer.box_size,
                    ctx.controls.ele_cutoff,
                ),
                &pr_reus,
            )?,
        )?;
    } else {
        for win in &windows {
            let win_dir = campaign_dir.join(format!("win_{}", win.id));
            fs::create_dir_all(&win_dir)?;
            write(&win_dir.join("system.gro"), &win.gro)?;
            if let Some(pdb) = &win.pdb {
                write(&win_dir.join("system.pdb"), pdb)?;
            }
            let restraint = ZRestraint {
                group1: &group1,
                group2: ctx.group2,
                function: func_str,
                direction: "Z2",
                constant: ctx.args.k,
                reference_angstrom: win.z_a,
            };
            write(
                &win_dir.join("min.inp"),
                &isolf::inp::minimization_with_restraint(
                    "../system.top",
                    "./system.gro",
                    "min",
                    ctx.min_params,
                    Some(&restraint),
                )?,
            )?;
            let eq_seed = rng.random_range(1..=MAX_ISEED);
            let pr_seed = rng.random_range(1..=MAX_ISEED);
            write(
                &win_dir.join("equil.inp"),
                &isolf::inp::dynamics_with_restraint(
                    ensemble,
                    "../system.top",
                    "./system.gro",
                    "./min.rst",
                    "equil",
                    &ctx.controls.equil.to_params(
                        ctx.controls.temperature,
                        eq_seed,
                        ctx.bilayer.box_size,
                        ctx.controls.ele_cutoff,
                    ),
                    Some(&restraint),
                )?,
            )?;
            write(
                &win_dir.join("prod.inp"),
                &isolf::inp::dynamics_with_restraint(
                    ensemble,
                    "../system.top",
                    "./system.gro",
                    "./equil.rst",
                    "prod",
                    &ctx.controls.prod.to_params(
                        ctx.controls.temperature,
                        pr_seed,
                        ctx.bilayer.box_size,
                        ctx.controls.ele_cutoff,
                    ),
                    Some(&restraint),
                )?,
            )?;
        }
    }

    let mbar_cfg = MbarConfig {
        cvfile: "./cv/{}.dat",
        fenefile: "./mbar.fene",
        pmffile: "./mbar.pmf",
        weightfile: Some("./weight_{}.dat"),
        temperature: ctx.controls.temperature,
        nblocks: ctx.args.mbar_blocks,
        grid_min_angstrom: ctx.grid_min_a,
        grid_max_angstrom: ctx.grid_max_a,
        grid_points: ctx.grid_points,
        constants: ctx.k_vec,
        references_angstrom: ctx.z_angstrom,
    };
    write(&campaign_dir.join("mbar.inp"), &isolf::inp::mbar(&mbar_cfg))?;

    let manifest = UmbrellaManifest {
        probe: probe_name,
        probe_kind: probe.kind,
        probe_atoms: probe.protein.beads.len(),
        mode: mode_label,
        ensemble: ctx.args.ensemble.label(),
        temperature_k: ctx.controls.temperature,
        k_kcal_mol_a2: ctx.args.k,
        restraint_function: func_str,
        restraint_direction: "Z2",
        group1: &group1,
        group2: ctx.group2,
        box_nm: ctx.bilayer.box_size,
        lipid_counts: &ctx.bilayer.lipid_counts,
        z_windows_nm: ctx.z_windows,
    };
    write(&campaign_dir.join("manifest.json"), &manifest.to_json())?;

    report::ok(
        probe_name,
        &format!(
            "{} windows (z = {:.2}..{:.2} nm, k = {} kcal/mol/A^2, {mode_label}) {} {} probe bead{} {} {} lipids",
            ctx.z_windows.len(),
            ctx.z_windows.first().copied().unwrap_or(0.0),
            ctx.z_windows.last().copied().unwrap_or(0.0),
            ctx.args.k,
            report::dot(),
            probe.protein.beads.len(),
            if probe.protein.beads.len() == 1 {
                ""
            } else {
                "s"
            },
            report::dot(),
            report::thousands(ctx.bilayer.total_lipids()),
        ),
    );
    Ok(())
}

/// `isolf umbrella setup`: build 1D bilayer umbrella-sampling windows (or REUS)
/// plus `lipid_protein.itp`, `mbar.inp`, and `manifest.json` for single amino-acid
/// beads or multi-residue PDB/mmCIF peptides/proteins.
pub(crate) fn run_umbrella_setup(args: UmbrellaSetupArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    let setup = validate_setup(&args, &force_field)?;

    let seed = args.general.seed.unwrap_or_else(|| rand::rng().random());
    let mut rng = StdRng::seed_from_u64(seed);
    let bilayer = load_or_build_bilayer(&args, &force_field, setup.min_box_z, &mut rng)?;

    let group2 = umbrella::bilayer_selection(&bilayer.lipid_counts);
    let z_angstrom: Vec<f64> = setup.z_windows.iter().map(|&z| z * 10.0).collect();
    let k_vec = vec![args.k; setup.z_windows.len()];
    let z_min_a = z_angstrom.iter().copied().fold(f64::INFINITY, f64::min);
    let z_max_a = z_angstrom.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let grid_min_a = z_min_a - args.mbar_grid_step;
    let grid_max_a = z_max_a + args.mbar_grid_step;
    let grid_points =
        (((grid_max_a - grid_min_a) / args.mbar_grid_step).round() as usize + 1).max(3);

    let isolf_itp = force_field.to_itp();
    let aicg2p_itp = render_aicg2p_itp();
    let min_params = setup.controls.min.to_params(
        bilayer.box_size,
        setup.controls.temperature,
        setup.controls.ele_cutoff,
    );
    let ctx = CampaignContext {
        args: &args,
        force_field: &force_field,
        bilayer: &bilayer,
        controls: &setup.controls,
        min_params: &min_params,
        z_windows: &setup.z_windows,
        z_angstrom: &z_angstrom,
        k_vec: &k_vec,
        group2: &group2,
        isolf_itp: &isolf_itp,
        aicg2p_itp: &aicg2p_itp,
        grid_min_a,
        grid_max_a,
        grid_points,
    };

    let multi_probe = setup.probes.len() > 1;
    for (probe, lp_params) in setup.probes.iter().zip(&setup.resolved_params) {
        let campaign_dir = if multi_probe {
            args.out.join(&probe.protein.name)
        } else {
            args.out.clone()
        };
        write_probe_campaign(&ctx, probe, lp_params, &campaign_dir, &mut rng)?;
    }

    if let Some(note) = ele_cutoff_note(bilayer.box_size, setup.controls.ele_cutoff) {
        report::warn(&note);
    }

    let out_dir = args.out.display().to_string();
    let cd = if out_dir == "." {
        String::new()
    } else {
        format!("cd {out_dir} && ")
    };
    let next = match args.mode {
        UmbrellaModeArg::Windows => format!(
            "run windows with atdyn (win_*/min.inp -> equil.inp -> prod.inp > prod.log), \
             then: isolf umbrella cv -o {out_dir} && {cd}mbar_analysis mbar.inp"
        ),
        UmbrellaModeArg::Reus => format!(
            "minimize min_*.inp, run mpirun -np {} atdyn equil_reus.inp -> prod_reus.inp, \
             then: isolf umbrella cv -o {out_dir} && {cd}mbar_analysis mbar.inp",
            setup.z_windows.len()
        ),
    };
    let prot_summary: Vec<(String, usize)> = setup
        .probes
        .iter()
        .map(|p| (p.protein.name.clone(), 1))
        .collect();
    let first_probe_beads = setup.probes.first().map_or(0, |p| p.protein.beads.len());
    report::summary(&report::Summary {
        identity: "umbrella",
        box_nm: bilayer.box_size,
        proteins: &prot_summary,
        lipids: &bilayer.lipid_counts,
        leaflets: None,
        beads: bilayer.particles.len() + first_probe_beads,
        next: Some(next),
    });
    Ok(())
}

/// `isolf umbrella cv`: manifest-free extraction of `RESTR_CVS001` from GENESIS
/// `atdyn` logs (and REUS `.rem` demuxing) into `cv/{}.dat` for `mbar_analysis`.
pub(crate) fn run_umbrella_cv(args: UmbrellaCvArgs) -> Result<(), Box<dyn Error>> {
    let targets = discover_cv_targets(&args.out, &args.log, &args.reus_prefix)?;
    report::header();
    report::section("Extract CV");

    for (dir, source) in targets {
        let label = dir
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("umbrella")
            .to_string();
        let cv_dir = dir.join("cv");
        fs::create_dir_all(&cv_dir)?;

        let (n_windows, total_samples) = match source {
            CvSource::Reus { count } => {
                let mut logs = Vec::with_capacity(count);
                let mut rems = Vec::with_capacity(count);
                for i in 1..=count {
                    let log_p = dir.join(format!("{}_{i}.log", args.reus_prefix));
                    let rem_p = dir.join(format!("{}_{i}.rem", args.reus_prefix));
                    logs.push(fs::read_to_string(&log_p).map_err(|e| {
                        format!("umbrella cv: failed to read {}: {e}", log_p.display())
                    })?);
                    rems.push(fs::read_to_string(&rem_p).map_err(|e| {
                        format!("umbrella cv: failed to read {}: {e}", rem_p.display())
                    })?);
                }
                let pairs: Vec<(&str, &str)> = logs
                    .iter()
                    .zip(&rems)
                    .map(|(l, r)| (l.as_str(), r.as_str()))
                    .collect();
                let demuxed = umbrella::demux_reus_cv(&pairs, args.discard)?;
                let mut samples = 0usize;
                for (idx, rows) in demuxed.iter().enumerate() {
                    samples += rows.len();
                    write(
                        &cv_dir.join(format!("{}.dat", idx + 1)),
                        &umbrella::render_cv_dat(rows),
                    )?;
                }
                (count, samples)
            }
            CvSource::Windows { count } => {
                let mut samples = 0usize;
                for i in 1..=count {
                    let log_p = dir.join(format!("win_{i}")).join(&args.log);
                    let log_txt = fs::read_to_string(&log_p).map_err(|e| {
                        format!("umbrella cv: failed to read {}: {e}", log_p.display())
                    })?;
                    let rows = umbrella::extract_cv_from_log(&log_txt, args.discard)?;
                    samples += rows.len();
                    write(
                        &cv_dir.join(format!("{i}.dat")),
                        &umbrella::render_cv_dat(&rows),
                    )?;
                }
                (count, samples)
            }
        };

        report::ok(
            &label,
            &format!(
                "{n_windows} windows {} {} samples -> {}",
                report::dot(),
                report::thousands(total_samples),
                cv_dir.display()
            ),
        );
    }
    Ok(())
}

fn write_params_file(
    target: &Path,
    label: &str,
    verb: &str,
    resolved: &LipidProteinParams,
) -> Result<(), Box<dyn Error>> {
    if let Some(parent) = target.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    write(target, &resolved.render_itp())?;
    report::ok(
        label,
        &format!(
            "{verb} {} ({} ISOLF, {} WCA, {} LJ)",
            target.display(),
            resolved.isolf.len(),
            resolved.wca.len(),
            resolved.lj.len()
        ),
    );
    Ok(())
}

/// `isolf params`: generate or update `lipid_protein.itp` for a standalone file,
/// an umbrella campaign, a multi-probe parent directory, or any membrane+protein
/// build directory.
pub(crate) fn run_params(args: ParamsArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    let base_params = match &args.params {
        Some(path) => {
            let content = fs::read_to_string(path)
                .map_err(|e| format!("--params: failed to read {}: {e}", path.display()))?;
            Some(LipidProteinParams::parse(&content)?)
        }
        None => None,
    };
    let overrides = ProbeParamOverrides {
        epsilon: args.eps,
        sigma: args.sigma,
        omega: args.omega,
        auto_head_wca: !args.no_auto_wca,
    };
    let out_is_file = args
        .out
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("itp"))
        || (args.out.exists() && args.out.is_file());

    report::header();
    report::section("Params");

    if let Some(probe_spec) = &args.probe {
        let probes = parse_umbrella_probes(probe_spec, None)?;
        let multi_subdirs = !out_is_file
            && probes.len() > 1
            && probes
                .iter()
                .all(|p| args.out.join(&p.protein.name).is_dir());
        if multi_subdirs {
            for probe in &probes {
                let target = args.out.join(&probe.protein.name).join("lipid_protein.itp");
                let resolved = LipidProteinParams::for_proteins(
                    &force_field,
                    &[&probe.protein],
                    base_params.as_ref(),
                    overrides,
                )?;
                write_params_file(&target, &probe.protein.name, "updated", &resolved)?;
            }
        } else {
            let target = if out_is_file {
                args.out.clone()
            } else {
                args.out.join("lipid_protein.itp")
            };
            let refs: Vec<_> = probes.iter().map(|p| &p.protein).collect();
            let resolved = LipidProteinParams::for_proteins(
                &force_field,
                &refs,
                base_params.as_ref(),
                overrides,
            )?;
            write_params_file(&target, "params", "wrote", &resolved)?;
        }
        return Ok(());
    }

    if out_is_file {
        let base = base_params.as_ref().ok_or(
            "isolf params: pass --probe <SPEC> or --params <FILE> when writing to a file path",
        )?;
        let aas: Vec<&str> = base.isolf.iter().map(|e| e.name.as_str()).collect();
        let resolved = LipidProteinParams::resolve(&force_field, &aas, Some(base), overrides)?;
        return write_params_file(&args.out, "params", "wrote", &resolved);
    }

    for DirResidues { dir, residues } in discover_protein_itp_dirs(&args.out)? {
        let label = dir.file_name().and_then(|s| s.to_str()).unwrap_or("params");
        let aa_refs: Vec<&str> = residues.iter().map(String::as_str).collect();
        let resolved =
            LipidProteinParams::resolve(&force_field, &aa_refs, base_params.as_ref(), overrides)?;
        write_params_file(&dir.join("lipid_protein.itp"), label, "updated", &resolved)?;
    }
    Ok(())
}

/// Find directories containing GENESIS umbrella/REUS logs (either `root` itself
/// or its immediate subdirectories).
fn discover_cv_targets(
    root: &Path,
    win_log: &str,
    reus_prefix: &str,
) -> Result<Vec<(PathBuf, CvSource)>, Box<dyn Error>> {
    if let Some(src) = umbrella::detect_cv_source(root, win_log, reus_prefix) {
        return Ok(vec![(root.to_path_buf(), src)]);
    }
    let mut found = Vec::new();
    if let Ok(entries) = fs::read_dir(root) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_dir()
                && let Some(src) = umbrella::detect_cv_source(&p, win_log, reus_prefix)
            {
                found.push((p, src));
            }
        }
    }
    found.sort_by(|a, b| a.0.cmp(&b.0));
    if found.is_empty() {
        return Err(format!(
            "no GENESIS umbrella logs found in {} (expected win_1/{win_log} or {reus_prefix}_1.log + .rem)",
            root.display()
        )
        .into());
    }
    Ok(found)
}

struct DirResidues {
    dir: PathBuf,
    residues: Vec<String>,
}

/// Discover directories containing protein `.itp` files and extract their
/// unique canonical amino-acid residue names (excluding `isolf.itp`,
/// `aicg2p.itp`, and `lipid_protein.itp`).
fn discover_protein_itp_dirs(root: &Path) -> Result<Vec<DirResidues>, Box<dyn Error>> {
    let direct = scan_dir_protein_residues(root);
    if !direct.is_empty() {
        return Ok(vec![DirResidues {
            dir: root.to_path_buf(),
            residues: direct,
        }]);
    }
    let mut found = Vec::new();
    if let Ok(entries) = fs::read_dir(root) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_dir() {
                let res = scan_dir_protein_residues(&p);
                if !res.is_empty() {
                    found.push(DirResidues {
                        dir: p,
                        residues: res,
                    });
                }
            }
        }
    }
    found.sort_by(|a, b| a.dir.cmp(&b.dir));
    if found.is_empty() {
        return Err(format!(
            "no protein .itp files found in {} (pass --probe <SPEC> to specify residues explicitly)",
            root.display()
        )
        .into());
    }
    Ok(found)
}

/// Extract unique canonical amino-acid names from the `[ atoms ]` section of
/// protein `.itp` files in `dir`.
fn scan_dir_protein_residues(dir: &Path) -> Vec<String> {
    let mut paths: Vec<PathBuf> = fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.extension().is_some_and(|e| e.eq_ignore_ascii_case("itp"))
                && !matches!(
                    p.file_name().and_then(|n| n.to_str()),
                    Some("isolf.itp" | "aicg2p.itp" | "lipid_protein.itp")
                )
        })
        .collect();
    paths.sort();

    let mut residues: Vec<String> = Vec::new();
    for path in paths {
        let Ok(content) = fs::read_to_string(&path) else {
            continue;
        };
        let mut in_atoms = false;
        for raw in content.lines() {
            let line = raw.split(';').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                let sec = line[1..line.len() - 1].trim();
                in_atoms = sec.eq_ignore_ascii_case("atoms");
                continue;
            }
            if !in_atoms {
                continue;
            }
            // GROMACS [ atoms ] columns: nr type resnr res atom cg charge mass
            if let Some(atom_type) = line.split_whitespace().nth(1)
                && let Some(canon) = canonical_aa(atom_type)
                && canon != "UNK"
                && !residues.iter().any(|r| r == canon)
            {
                residues.push(canon.to_string());
            }
        }
    }
    residues
}
