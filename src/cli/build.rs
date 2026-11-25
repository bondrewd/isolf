//! Lipid geometry build runners (`membrane`, `vesicle`, `tubule`) and shared
//! output writers (`.gro`, `.top`, `.inp`, `.vmd`, `.gif`).

use std::error::Error;
use std::f64::consts::TAU;
use std::fs;
use std::path::Path;

use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;

use isolf::composition::Composition;
use isolf::embed::{
    EmbeddedSystem, embed_membrane, embed_membrane_recorded, embed_tubule, embed_tubule_recorded,
    embed_vesicle, embed_vesicle_recorded,
};
use isolf::error::BuildError;
use isolf::force_field::{ForceField, LipidProteinParams, ProbeParamOverrides};
use isolf::membrane::{BuildOptions, DEFAULT_NAME, Particle, RelaxFrame, Sizing};
use isolf::protein::params::render_aicg2p_itp;
use isolf::tubule::TubuleRadius;
use isolf::vesicle::VesicleRadius;
use isolf::vmd::{Layout, VmdOptions, VmdSource};

use crate::{anim, report, runlog};

use super::args::{
    CenterArg, CommonBuildArgs, GifModeArg, LipidColorScope, MembraneArgs, PhaseArg, TubuleArgs,
    VesicleArgs, VmdArgs,
};
use super::lipids::enrich_build_error;
use super::parse::{
    ControlSpecs, parse_box_spec, parse_count_spec, parse_protein_specs, parse_radius_spec,
    resolve_compositions, validate_padding,
};
use super::{MAX_ISEED, reproduce_command, write};

/// One optional structure-file format: whether it was requested, its extension,
/// and how to render it from the built system.
type StructureFormat = (bool, &'static str, fn(&EmbeddedSystem) -> String);

/// What is being built, which selects the production control file, the camera,
/// and the default output naming.
#[derive(Clone, Copy, PartialEq)]
pub(crate) enum Geometry {
    Membrane,
    Vesicle,
    Tubule,
}

impl Geometry {
    fn label(self) -> &'static str {
        match self {
            Geometry::Membrane => "membrane",
            Geometry::Vesicle => "vesicle",
            Geometry::Tubule => "tubule",
        }
    }

    fn title(self) -> &'static str {
        match self {
            Geometry::Membrane => DEFAULT_NAME,
            Geometry::Vesicle => "CG vesicle model",
            Geometry::Tubule => "CG tubule model",
        }
    }

    fn layout(self) -> Layout {
        match self {
            Geometry::Membrane => Layout::Membrane,
            Geometry::Vesicle => Layout::Vesicle,
            Geometry::Tubule => Layout::Tubule,
        }
    }

    fn is_npt(self) -> bool {
        self == Geometry::Membrane
    }
}

/// Emit a transparent warning/note whenever `cg_cutoffdist_ele` is auto-scaled
/// below GENESIS's built-in 52.0 Å default to satisfy `box > 2 * cutoff`.
pub(crate) fn ele_cutoff_note(box_size: [f64; 3], explicit: Option<f64>) -> Option<String> {
    if explicit.is_some() {
        return None;
    }
    let (rc, rp) = isolf::inp::ele_cutoff_for_box(box_size)?;
    let min_nm = box_size[0].min(box_size[1]).min(box_size[2]);
    Some(format!(
        "small box (min axis {min_nm:.2} nm < 12.0 nm): auto-set cg_cutoffdist_ele = {rc:.1} A, \
         cg_pairlistdist_ele = {rp:.1} A for GENESIS PBC (override with --ele-cutoff)"
    ))
}

fn load_lipid_protein_opt(
    path: Option<&Path>,
) -> Result<Option<LipidProteinParams>, Box<dyn Error>> {
    match path {
        Some(p) => {
            let content = fs::read_to_string(p)
                .map_err(|e| format!("--lipid-protein: failed to read {}: {e}", p.display()))?;
            Ok(Some(LipidProteinParams::parse(&content)?))
        }
        None => Ok(None),
    }
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

pub(crate) fn run_membrane(args: MembraneArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    validate_padding(args.padding)?;

    let (upper, lower) = resolve_compositions(
        &args.bilayer.lipids,
        &args.bilayer.upper,
        &args.bilayer.lower,
        "--upper",
        "--lower",
    )?;
    let counts = parse_count_spec(&args.bilayer.count)?;
    let box_xy = parse_box_spec(&args.bilayer.box_size)?;
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
    let lp_params = load_lipid_protein_opt(args.common.lipid_protein.as_deref())?;
    let (controls, temp_note) = args
        .common
        .controls
        .prepare(args.common.want_inp(), Geometry::Membrane.is_npt())?;

    let seed = args
        .common
        .general
        .seed
        .unwrap_or_else(|| rand::rng().random());
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
        lp_params.as_ref(),
        &args.common,
        Geometry::Membrane,
        args.center,
        &controls,
        &notes,
        seed,
        &mut rng,
    )
}

pub(crate) fn run_vesicle(args: VesicleArgs) -> Result<(), Box<dyn Error>> {
    let force_field = ForceField::isolf();
    validate_padding(args.padding)?;

    let (outer, inner) =
        resolve_compositions(&args.lipids, &args.outer, &args.inner, "--outer", "--inner")?;
    let radius = parse_radius_spec(&args.radius)?
        .ok_or("missing vesicle radius: pass --radius <R> (e.g. 20, ri=15, or ri=14,ro=20)")?;
    isolf::membrane::area_per_lipid(&force_field, &outer)
        .map_err(|e| enrich_build_error(&force_field, e))?;
    isolf::membrane::area_per_lipid(&force_field, &inner)
        .map_err(|e| enrich_build_error(&force_field, e))?;
    let protein_specs = parse_protein_specs(&args.common.proteins)?;
    let lp_params = load_lipid_protein_opt(args.common.lipid_protein.as_deref())?;
    let (controls, temp_note) = args
        .common
        .controls
        .prepare(args.common.want_inp(), Geometry::Vesicle.is_npt())?;

    let seed = args
        .common
        .general
        .seed
        .unwrap_or_else(|| rand::rng().random());
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
        lp_params.as_ref(),
        &args.common,
        Geometry::Vesicle,
        args.center,
        &controls,
        &notes,
        seed,
        &mut rng,
    )
}

pub(crate) fn run_tubule(args: TubuleArgs) -> Result<(), Box<dyn Error>> {
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
    let lp_params = load_lipid_protein_opt(args.common.lipid_protein.as_deref())?;
    let (controls, temp_note) = args
        .common
        .controls
        .prepare(args.common.want_inp(), Geometry::Tubule.is_npt())?;

    let seed = args
        .common
        .general
        .seed
        .unwrap_or_else(|| rand::rng().random());
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
        lp_params.as_ref(),
        &args.common,
        Geometry::Tubule,
        args.center,
        &controls,
        &notes,
        seed,
        &mut rng,
    )
}

#[allow(clippy::too_many_arguments)]
fn finish_build(
    mut system: EmbeddedSystem,
    gif_frames: Vec<RelaxFrame>,
    gif_bounds: [f64; 4],
    force_field: &ForceField,
    lp_base: Option<&LipidProteinParams>,
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
    if lp_base.is_some() && !system.has_proteins() {
        report::warn("--lipid-protein has no effect without embedded proteins (-p)");
    }
    if common.want_inp()
        && let Some(note) = ele_cutoff_note(system.box_size, controls.ele_cutoff)
    {
        report::warn(&note);
    }

    let base = common
        .name
        .clone()
        .unwrap_or_else(|| geometry.label().to_string());

    write_outputs(
        &mut system,
        force_field,
        lp_base,
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
    let command = reproduce_command(seed, common.general.seed.is_some());
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

/// Resolve the flat-membrane sizing, matching the two leaflets' areas in count
/// mode unless `no_balance` is set. Also validates all lipid names up front.
pub(crate) fn resolve_membrane_sizing(
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

/// Per-leaflet lipid composition (excluding protein particles).
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
    lp_base: Option<&LipidProteinParams>,
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
        (common.formats.pdb, "pdb", EmbeddedSystem::to_pdb),
        (
            common.formats.psf || want_vmd,
            "psf",
            EmbeddedSystem::to_psf,
        ),
        (common.formats.crd, "crd", EmbeddedSystem::to_crd),
        (common.formats.cif, "cif", EmbeddedSystem::to_cif),
    ];
    for (enabled, ext, render) in formats {
        if enabled {
            write(&named(ext), &render(system))?;
        }
    }

    let include_lp = lp_base.is_some() && system.has_proteins();
    if common.top || want_inp || include_lp {
        write(&named("top"), &system.to_top_with_lipid_protein(include_lp))?;
        write(&common.out.join("isolf.itp"), &force_field.to_itp())?;
        if system.has_proteins() {
            write(&common.out.join("aicg2p.itp"), &render_aicg2p_itp())?;
            if let Some(lp) = lp_base {
                let refs: Vec<_> = system.protein_templates.iter().collect();
                let resolved = LipidProteinParams::for_proteins(
                    force_field,
                    &refs,
                    Some(lp),
                    ProbeParamOverrides {
                        auto_head_wca: true,
                        ..ProbeParamOverrides::default()
                    },
                )?;
                write(
                    &common.out.join("lipid_protein.itp"),
                    &resolved.render_itp(),
                )?;
            }
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
        write_control_files(
            &common.out,
            base,
            geometry.is_npt(),
            box_size,
            controls,
            rng,
        )?;
    }
    if want_vmd {
        write_vmd_scripts(
            &common.out,
            base,
            geometry.layout(),
            geometry.is_npt(),
            want_inp,
            matches!(center, CenterArg::Origin),
            &system.particles,
            force_field,
            &common.vmd_opts,
            common.vmd_color_by,
            controls,
        )?;
    }
    Ok(())
}

/// Write the GENESIS control files: minimize, then equilibrate and run
/// production under the requested ensemble (`is_npt`: NPT for flat membranes,
/// NVT for vesicles, tubules, and standalone proteins).
pub(crate) fn write_control_files(
    out: &Path,
    base: &str,
    is_npt: bool,
    box_size: [f64; 3],
    controls: &ControlSpecs,
    rng: &mut StdRng,
) -> Result<(), Box<dyn Error>> {
    let min_params = controls
        .min
        .to_params(box_size, controls.temperature, controls.ele_cutoff);
    write(
        &out.join("min.inp"),
        &isolf::inp::minimization(base, &min_params)?,
    )?;

    let prod_seed = rng.random_range(1..=MAX_ISEED);
    let equil_seed = rng.random_range(1..=MAX_ISEED);
    let production = controls.prod.to_params(
        controls.temperature,
        prod_seed,
        box_size,
        controls.ele_cutoff,
    );
    let equilibration = controls.equil.to_params(
        controls.temperature,
        equil_seed,
        box_size,
        controls.ele_cutoff,
    );

    if is_npt {
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
/// when `want_inp` is enabled).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_vmd_scripts(
    out: &Path,
    base: &str,
    layout: Layout,
    is_npt: bool,
    want_inp: bool,
    recenter_to_box: bool,
    particles: &[Particle],
    force_field: &ForceField,
    vmd: &VmdArgs,
    color_by: Option<LipidColorScope>,
    controls: &ControlSpecs,
) -> Result<(), Box<dyn Error>> {
    let background = vmd.vmd_background.as_deref().unwrap_or("white");
    let material = vmd.vmd_material.as_deref().unwrap_or("AOChalky");
    let options = VmdOptions {
        surface: vmd.vmd_surface,
        cutaway: vmd.vmd_cutaway,
        render: vmd.vmd_render,
        background,
        material,
        smooth: vmd.vmd_smooth.unwrap_or(5),
        lipid_coloring: match color_by.unwrap_or(LipidColorScope::Head) {
            LipidColorScope::Head => isolf::vmd::LipidColoring::Head,
            LipidColorScope::Bead => isolf::vmd::LipidColoring::Role,
            LipidColorScope::Whole => isolf::vmd::LipidColoring::Whole,
        },
        recenter_to_box,
    };
    let script =
        |source| isolf::vmd::vmd_script(force_field, particles, base, layout, source, &options);

    write(
        &out.join(format!("{base}.vmd")),
        &script(VmdSource::Structure),
    )?;

    if want_inp {
        let equil_name = if is_npt { "npt" } else { "nvt" };
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
}
