use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn isolf_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_isolf"))
}

fn temp_out(label: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("isolf_cli_{label}_{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    dir
}

/// Write a synthetic 20-residue OPM-aligned transmembrane α-helix PDB centered
/// at `z = 0` (spans z ≈ -1.4 to +1.4 nm) for self-contained CLI tests.
fn write_sample_helix_pdb(path: &std::path::Path) {
    let mut s = String::new();
    let mut serial = 1;
    let seq = [
        "ALA", "LEU", "VAL", "ILE", "LEU", "PHE", "ALA", "TRP", "LEU", "LYS", "GLU", "LEU", "VAL",
        "ALA", "ILE", "PHE", "LEU", "ARG", "ASP", "GLY",
    ];
    for (i, &res) in seq.iter().enumerate() {
        let idx = i + 1;
        let theta = i as f64 * (100.0_f64).to_radians();
        let r = 2.3;
        let z = (i as f64 - 9.5) * 1.5;
        let ca = [r * theta.cos(), r * theta.sin(), z];
        let n = [r * (theta - 0.3).cos(), r * (theta - 0.3).sin(), z - 0.6];
        let c = [r * (theta + 0.3).cos(), r * (theta + 0.3).sin(), z + 0.6];
        let o = [
            (r + 0.4) * (theta + 0.3).cos(),
            (r + 0.4) * (theta + 0.3).sin(),
            z + 1.3,
        ];
        for (aname, p) in [("N", n), ("CA", ca), ("C", c), ("O", o)] {
            let _ = writeln!(
                s,
                "ATOM  {:>5} {:<4} {:<3} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00",
                serial, aname, res, idx, p[0], p[1], p[2]
            );
            serial += 1;
        }
    }
    s.push_str("TER\nEND\n");
    fs::write(path, s).unwrap();
}

/// Write the same helix as an mmCIF file to verify mmCIF parsing end-to-end.
fn write_sample_helix_cif(path: &std::path::Path) {
    let mut s = String::from(
        "data_helix\n#\nloop_\n\
         _atom_site.group_PDB\n\
         _atom_site.id\n\
         _atom_site.label_atom_id\n\
         _atom_site.label_alt_id\n\
         _atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n\
         _atom_site.label_seq_id\n\
         _atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n\
         _atom_site.Cartn_z\n\
         _atom_site.pdbx_PDB_model_num\n",
    );
    let mut serial = 1;
    for i in 0..12 {
        let idx = i + 1;
        let theta = i as f64 * (100.0_f64).to_radians();
        let r = 2.3;
        let z = (i as f64 - 5.5) * 1.5;
        let ca = [r * theta.cos(), r * theta.sin(), z];
        let n = [r * (theta - 0.3).cos(), r * (theta - 0.3).sin(), z - 0.6];
        let c = [r * (theta + 0.3).cos(), r * (theta + 0.3).sin(), z + 0.6];
        let o = [
            (r + 0.4) * (theta + 0.3).cos(),
            (r + 0.4) * (theta + 0.3).sin(),
            z + 1.3,
        ];
        for (aname, p) in [("N", n), ("CA", ca), ("C", c), ("O", o)] {
            let _ = writeln!(
                s,
                "ATOM {} {} . LEU A {} {:.3} {:.3} {:.3} 1",
                serial, aname, idx, p[0], p[1], p[2]
            );
            serial += 1;
        }
    }
    s.push_str("#\n");
    fs::write(path, s).unwrap();
}

#[test]
fn membrane_phase_raft_writes_full_genesis_vmd_and_gif_bundle() {
    let out = temp_out("membrane_raft");
    let status = isolf_bin()
        .args([
            "membrane",
            "-L",
            "DPPC=1,DSPC=1;POPC=2,POPS=1",
            "-n",
            "80",
            "--phase",
            "raft",
            "--anti-register",
            "--inp",
            "--vmd",
            "--gif",
            "--seed",
            "7",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&out)
        .output()
        .expect("run isolf membrane");

    assert!(
        status.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    for file in [
        "membrane.gro",
        "membrane.log",
        "membrane.top",
        "isolf.itp",
        "membrane.psf",
        "membrane.vmd",
        "membrane.gif",
        "min.inp",
        "npt.inp",
        "pro.inp",
        "min.vmd",
        "npt.vmd",
        "pro.vmd",
    ] {
        let p = out.join(file);
        assert!(p.exists(), "missing expected output {file}");
        assert!(fs::metadata(&p).unwrap().len() > 0, "empty output {file}");
    }
    let top = fs::read_to_string(out.join("membrane.top")).unwrap();
    for lipid in ["DPPC", "DSPC", "POPC", "POPS"] {
        assert!(top.contains(lipid), "topology missing {lipid}");
    }
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn vesicle_phase_stripe_uses_nvt_controls() {
    let out = temp_out("vesicle_stripe");
    let output = isolf_bin()
        .args([
            "vesicle",
            "-L",
            "DPPC=1,POPC=1",
            "-r",
            "8",
            "--phase",
            "stripe",
            "--inp",
            "--vmd",
            "--seed",
            "3",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&out)
        .output()
        .expect("run isolf vesicle");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(out.join("vesicle.gro").exists());
    assert!(out.join("nvt.inp").exists());
    assert!(out.join("pro.inp").exists());
    assert!(out.join("nvt.vmd").exists());
    assert!(!out.join("npt.inp").exists());
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn tubule_phase_raft_writes_periodic_z_gro_nvt_vmd_and_gif() {
    let out = temp_out("tubule_raft");
    let output = isolf_bin()
        .args([
            "tubule",
            "-L",
            "DPPC=1,POPC=3",
            "-r",
            "8",
            "-l",
            "16",
            "--phase",
            "raft",
            "--anti-register",
            "--inp",
            "--vmd",
            "--gif",
            "--seed",
            "5",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&out)
        .output()
        .expect("run isolf tubule");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    for file in [
        "tubule.gro",
        "tubule.log",
        "tubule.top",
        "isolf.itp",
        "tubule.psf",
        "tubule.vmd",
        "tubule.gif",
        "min.inp",
        "nvt.inp",
        "pro.inp",
    ] {
        let p = out.join(file);
        assert!(p.exists(), "missing expected output {file}");
        assert!(fs::metadata(&p).unwrap().len() > 0, "empty output {file}");
    }
    assert!(!out.join("npt.inp").exists());

    let gro = fs::read_to_string(out.join("tubule.gro")).unwrap();
    let box_line = gro.lines().rfind(|l| !l.trim().is_empty()).unwrap();
    let dims: Vec<f64> = box_line
        .split_whitespace()
        .map(|s| s.parse::<f64>().unwrap())
        .collect();
    assert_eq!(dims.len(), 3);
    // 2 * R_outer (16) + 2 * default padding (20) = 36 nm in x/y; periodic z = 16 nm.
    assert!((dims[0] - 36.0).abs() < 1e-3, "box x={}", dims[0]);
    assert!((dims[1] - 36.0).abs() < 1e-3, "box y={}", dims[1]);
    assert!((dims[2] - 16.0).abs() < 1e-3, "box z={}", dims[2]);

    let vmd = fs::read_to_string(out.join("tubule.vmd")).unwrap();
    assert!(vmd.contains("tubule"), "tubule.vmd should identify layout");
    assert!(vmd.contains("display projection perspective"));
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn lipids_subcommand_and_diagnostics_work() {
    let list = isolf_bin()
        .args(["lipids", "--no-color", "--ascii"])
        .output()
        .expect("run isolf lipids");
    assert!(list.status.success());
    let stdout = String::from_utf8_lossy(&list.stdout);
    assert!(stdout.contains("35 species"));
    assert!(stdout.contains("POPC"));
    assert!(stdout.contains("DLPA"));
    assert!(stdout.contains("SOPS"));

    // Unknown lipid gives edit-distance suggestions.
    let bad = isolf_bin()
        .args(["membrane", "-L", "POPX", "-n", "16", "--no-color"])
        .output()
        .expect("run isolf with typo");
    assert!(!bad.status.success());
    let stderr = String::from_utf8_lossy(&bad.stderr);
    assert!(stderr.contains("POPC"), "stderr: {stderr}");

    // Too-small tubule radius fails cleanly.
    let small = isolf_bin()
        .args(["tubule", "-L", "POPC", "-r", "2", "-l", "15", "--no-color"])
        .output()
        .expect("run isolf with small tubule");
    assert!(!small.status.success());
    let stderr = String::from_utf8_lossy(&small.stderr);
    assert!(stderr.contains("bilayer thickness"), "stderr: {stderr}");

    // --anti-register with --phase random emits a warning rather than silently doing nothing.
    let out = temp_out("warn_antireg");
    let warn = isolf_bin()
        .args([
            "membrane",
            "-L",
            "DPPC=1,POPC=1",
            "-n",
            "32",
            "--anti-register",
            "--seed",
            "1",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&out)
        .output()
        .expect("run isolf anti-register warning");
    assert!(warn.status.success());
    let combined = format!(
        "{}{}",
        String::from_utf8_lossy(&warn.stdout),
        String::from_utf8_lossy(&warn.stderr)
    );
    assert!(
        combined.contains("--anti-register"),
        "expected warning in output: {combined}"
    );
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn protein_subcommand_builds_aicg2p_from_pdb_and_cif() {
    let out = temp_out("protein_standalone");
    fs::create_dir_all(&out).unwrap();
    let pdb_path = out.join("helix.pdb");
    let cif_path = out.join("helix.cif");
    let charge_path = out.join("charges.dat");
    write_sample_helix_pdb(&pdb_path);
    write_sample_helix_cif(&cif_path);
    fs::write(&charge_path, "# RESPAC override\n1 0.5\n").unwrap();

    let run_dir = out.join("pdb_out");
    let output = isolf_bin()
        .args(["protein"])
        .arg(&pdb_path)
        .args(["--idr-local", "1-3", "--idr-nonlocal", "1-2", "--charges"])
        .arg(&charge_path)
        .args([
            "--inp",
            "--vmd",
            "--pdb",
            "--crd",
            "--cif",
            "--seed",
            "7",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&run_dir)
        .output()
        .expect("run isolf protein on pdb");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    for file in [
        "helix_cg.gro",
        "helix_cg.itp",
        "aicg2p.itp",
        "helix_cg.top",
        "helix_cg.log",
        "helix_cg.pdb",
        "helix_cg.psf",
        "helix_cg.crd",
        "helix_cg.cif",
        "helix_cg.vmd",
        "min.inp",
        "nvt.inp",
        "pro.inp",
    ] {
        let p = run_dir.join(file);
        assert!(p.exists(), "missing expected protein output {file}");
        assert!(fs::metadata(&p).unwrap().len() > 0, "empty output {file}");
    }
    let itp = fs::read_to_string(run_dir.join("helix_cg.itp")).unwrap();
    assert!(itp.contains("0.500"), "charge override missing in itp");
    assert!(itp.contains("[ pairs ]"));
    let vmd = fs::read_to_string(run_dir.join("helix_cg.vmd")).unwrap();
    assert!(vmd.contains("Tube 1.8 12.0"));
    assert!(vmd.contains("name CA"));

    // mmCIF input with Clementi Gō model.
    let cif_dir = out.join("cif_out");
    let cif_out = isolf_bin()
        .args(["protein"])
        .arg(&cif_path)
        .args([
            "--model",
            "clementi",
            "--center",
            "origin",
            "--no-color",
            "-o",
        ])
        .arg(&cif_dir)
        .output()
        .expect("run isolf protein on cif");
    assert!(
        cif_out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&cif_out.stderr)
    );
    let cif_itp = fs::read_to_string(cif_dir.join("helix_cg.itp")).unwrap();
    assert!(cif_itp.contains("8.3680E+04"), "expected Clementi bond K");
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn embeds_proteins_in_membrane_vesicle_and_tubule_via_cli() {
    let out = temp_out("protein_embed");
    fs::create_dir_all(&out).unwrap();
    let pdb_path = out.join("tm.pdb");
    write_sample_helix_pdb(&pdb_path);

    // 1. Flat membrane with 2 embedded copies + --inp + --vmd.
    let mem_dir = out.join("mem");
    let prot_spec = format!("{},copies=2,clearance=0.2", pdb_path.display());
    let mem = isolf_bin()
        .args([
            "membrane",
            "-L",
            "POPC=3,POPS=1",
            "-n",
            "80",
            "-p",
            &prot_spec,
            "--inp",
            "--vmd",
            "--seed",
            "11",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&mem_dir)
        .output()
        .expect("run isolf membrane -p");
    assert!(
        mem.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&mem.stderr)
    );
    for file in [
        "membrane.gro",
        "membrane.top",
        "membrane.psf",
        "membrane.vmd",
        "isolf.itp",
        "aicg2p.itp",
        "tm_cg.itp",
        "min.inp",
        "npt.inp",
        "pro.inp",
    ] {
        assert!(
            mem_dir.join(file).exists(),
            "missing membrane output {file}"
        );
    }
    let top = fs::read_to_string(mem_dir.join("membrane.top")).unwrap();
    assert!(top.contains("#include \"./aicg2p.itp\""));
    assert!(top.contains("#include \"./isolf.itp\""));
    assert!(top.contains("#include \"./tm_cg.itp\""));
    assert!(top.contains("tm_cg 2"));
    let vmd = fs::read_to_string(mem_dir.join("membrane.vmd")).unwrap();
    assert!(vmd.contains("Tube 1.8 12.0"));
    assert!(vmd.contains("name PHO"));

    // 2. Vesicle with embedded protein.
    let ves_dir = out.join("ves");
    let ves_spec = format!("{},z=center,copies=2", pdb_path.display());
    let ves = isolf_bin()
        .args([
            "vesicle",
            "-L",
            "POPC",
            "-r",
            "8",
            "-p",
            &ves_spec,
            "--top",
            "--vmd",
            "--seed",
            "4",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&ves_dir)
        .output()
        .expect("run isolf vesicle -p");
    assert!(
        ves.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ves.stderr)
    );
    assert!(ves_dir.join("aicg2p.itp").exists());
    assert!(ves_dir.join("tm_cg.itp").exists());

    // 3. Tubule with embedded protein.
    let tub_dir = out.join("tub");
    let tub = isolf_bin()
        .args([
            "tubule",
            "-L",
            "POPC",
            "-r",
            "8",
            "-l",
            "16",
            "-p",
            &ves_spec,
            "--top",
            "--vmd",
            "--seed",
            "4",
            "--no-color",
            "--ascii",
            "-o",
        ])
        .arg(&tub_dir)
        .output()
        .expect("run isolf tubule -p");
    assert!(
        tub.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&tub.stderr)
    );
    assert!(tub_dir.join("aicg2p.itp").exists());
    assert!(tub_dir.join("tm_cg.itp").exists());

    let _ = fs::remove_dir_all(&out);
}
