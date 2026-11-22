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
