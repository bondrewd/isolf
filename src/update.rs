//! Self-management: update the running binary, or uninstall it.
//!
//! The `update` subcommand resolves the latest release of `bondrewd/isolf` and, if it is newer
//! than this build, downloads the archive for the current platform, verifies its
//! sha256, and swaps it in over the running executable. macOS uses the universal
//! binary, Linux the static musl build, and Windows the msvc build (a `.zip`,
//! swapped in with the move-aside dance a running `.exe` needs). The archive
//! layout and asset names match `.github/workflows/release.yml`.
//!
//! The `uninstall` subcommand deletes the running binary: Unix unlinks it directly, Windows
//! moves the locked `.exe` aside and schedules its deletion once the process exits.

use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const REPO: &str = "bondrewd/isolf";
const USER_AGENT: &str = concat!("isolf/", env!("CARGO_PKG_VERSION"));

/// Update the binary in place. The entry point for the `update` subcommand.
pub fn run() -> Result<(), Box<dyn Error>> {
    let current = env!("CARGO_PKG_VERSION");
    let latest = latest_version()?;
    if !is_newer(&latest, current) {
        println!("isolf {current} is already up to date");
        return Ok(());
    }
    let target = target_triple().ok_or_else(|| -> Box<dyn Error> {
        format!(
            "automatic update is not supported on this platform; \
             download the latest from https://github.com/{REPO}/releases"
        )
        .into()
    })?;
    println!("updating isolf {current} -> {latest}");

    // Windows ships a .zip holding isolf.exe; the Unix targets a .tar.gz / isolf.
    let (ext, bin_name) = if cfg!(windows) {
        ("zip", "isolf.exe")
    } else {
        ("tar.gz", "isolf")
    };
    let tag = format!("v{latest}");
    let base = format!("isolf-{tag}-{target}");
    let release = format!("https://github.com/{REPO}/releases/download/{tag}");

    // Stage alongside the target so the final swap stays on one filesystem.
    let exe = std::env::current_exe()?;
    let dir = exe.parent().ok_or("cannot locate the install directory")?;
    let archive = dir.join(format!(".isolf-update.{ext}"));
    let stage = dir.join(".isolf-update.d");
    let _cleanup = CleanUp(vec![archive.clone(), stage.clone()]);

    println!("downloading {base}.{ext}");
    download(&format!("{release}/{base}.{ext}"), &archive)?;
    verify(&archive, &format!("{release}/{base}.sha256"))?;

    let _ = fs::remove_dir_all(&stage);
    fs::create_dir_all(&stage)?;
    // `tar` unpacks the .tar.gz (auto-detecting gzip) and, via the bsdtar bundled
    // with Windows 10+, the .zip.
    let extract = if cfg!(windows) { "-xf" } else { "-xzf" };
    let status = Command::new("tar")
        .arg(extract)
        .arg(&archive)
        .arg("-C")
        .arg(&stage)
        .status()?;
    if !status.success() {
        return Err("failed to extract the downloaded archive".into());
    }
    let new_bin = stage.join(bin_name);
    if !new_bin.exists() {
        return Err("the downloaded archive did not contain an isolf binary".into());
    }
    set_executable(&new_bin)?;
    replace_exe(&new_bin, &exe).map_err(|e| -> Box<dyn Error> {
        format!(
            "could not replace {} ({e}); you may lack write permission there",
            exe.display()
        )
        .into()
    })?;
    println!("updated to isolf {latest}");
    Ok(())
}

/// Remove the running binary from the system. The entry point for the `uninstall` subcommand.
/// isolf keeps no global config or cache, so deleting the executable is the whole
/// of it.
pub fn uninstall() -> Result<(), Box<dyn Error>> {
    let exe = std::env::current_exe()?;
    remove_exe(&exe).map_err(|e| -> Box<dyn Error> {
        format!(
            "could not remove {} ({e}); you may lack write permission there",
            exe.display()
        )
        .into()
    })?;
    println!("uninstalled isolf ({})", exe.display());
    Ok(())
}

/// The release asset target triple for the platform this binary was built for.
fn target_triple() -> Option<&'static str> {
    if cfg!(target_os = "macos") {
        Some("universal-apple-darwin")
    } else if cfg!(target_os = "linux") && cfg!(target_arch = "x86_64") {
        Some("x86_64-unknown-linux-musl")
    } else if cfg!(target_os = "linux") && cfg!(target_arch = "aarch64") {
        Some("aarch64-unknown-linux-musl")
    } else if cfg!(target_os = "windows") {
        // The only Windows build; it runs on arm64 Windows through emulation.
        Some("x86_64-pc-windows-msvc")
    } else {
        None
    }
}

/// Resolve the latest release version (without the leading `v`) from GitHub.
fn latest_version() -> Result<String, Box<dyn Error>> {
    let url = format!("https://api.github.com/repos/{REPO}/releases/latest");
    let mut resp = ureq::get(&url).header("User-Agent", USER_AGENT).call()?;
    let body = resp.body_mut().read_to_string()?;
    let tag =
        json_field(&body, "tag_name").ok_or("could not read the latest version from GitHub")?;
    Ok(tag.trim_start_matches('v').to_string())
}

/// Download `url` to `path`.
fn download(url: &str, path: &Path) -> Result<(), Box<dyn Error>> {
    let mut resp = ureq::get(url).header("User-Agent", USER_AGENT).call()?;
    let bytes = resp
        .body_mut()
        .with_config()
        .limit(256 * 1024 * 1024)
        .read_to_vec()?;
    fs::write(path, bytes)?;
    Ok(())
}

/// Verify `archive` against the sha256 published at `sum_url`. Best-effort: with
/// no sha256 tool on the system it warns and proceeds (the download was HTTPS).
fn verify(archive: &Path, sum_url: &str) -> Result<(), Box<dyn Error>> {
    let mut resp = ureq::get(sum_url).header("User-Agent", USER_AGENT).call()?;
    let published = resp.body_mut().read_to_string()?;
    let expected = published
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_lowercase();
    match sha256(archive) {
        Some(actual) if actual == expected => Ok(()),
        Some(actual) => Err(format!("checksum mismatch: expected {expected}, got {actual}").into()),
        None => {
            eprintln!("note: no sha256 tool found, skipping checksum verification");
            Ok(())
        }
    }
}

/// The sha256 of `path` via a system tool, or `None` if none is available.
fn sha256(path: &Path) -> Option<String> {
    let tools: [(&str, &[&str]); 2] = [("sha256sum", &[]), ("shasum", &["-a", "256"])];
    for (cmd, args) in tools {
        let Ok(out) = Command::new(cmd).args(args).arg(path).output() else {
            continue;
        };
        if !out.status.success() {
            continue;
        }
        let text = String::from_utf8_lossy(&out.stdout);
        if let Some(hash) = text.split_whitespace().next() {
            return Some(hash.to_lowercase());
        }
    }
    // Windows: `certutil -hashfile <path> SHA256` prints the hash on its own line.
    if let Ok(out) = Command::new("certutil")
        .arg("-hashfile")
        .arg(path)
        .arg("SHA256")
        .output()
        && out.status.success()
    {
        for line in String::from_utf8_lossy(&out.stdout).lines() {
            let hex: String = line.split_whitespace().collect::<String>().to_lowercase();
            if hex.len() == 64 && hex.bytes().all(|b| b.is_ascii_hexdigit()) {
                return Some(hex);
            }
        }
    }
    None
}

/// Extract a string field from a flat JSON object by key.
fn json_field(body: &str, key: &str) -> Option<String> {
    let needle = format!("\"{key}\"");
    let after = &body[body.find(&needle)? + needle.len()..];
    let after = &after[after.find(':')? + 1..];
    let start = after.find('"')? + 1;
    let end = after[start..].find('"')?;
    Some(after[start..start + end].to_string())
}

/// Is `latest` a newer semantic version than `current`?
fn is_newer(latest: &str, current: &str) -> bool {
    semver(latest) > semver(current)
}

fn semver(v: &str) -> (u64, u64, u64) {
    let mut p = v
        .trim_start_matches('v')
        .split(['.', '-', '+'])
        .map(|s| s.parse().unwrap_or(0));
    (
        p.next().unwrap_or(0),
        p.next().unwrap_or(0),
        p.next().unwrap_or(0),
    )
}

#[cfg(unix)]
fn set_executable(path: &Path) -> std::io::Result<()> {
    use std::os::unix::fs::PermissionsExt;
    let mut perms = fs::metadata(path)?.permissions();
    perms.set_mode(0o755);
    fs::set_permissions(path, perms)
}

#[cfg(not(unix))]
fn set_executable(_path: &Path) -> std::io::Result<()> {
    Ok(())
}

/// Swap the new binary in over the running executable. On Unix a rename over the
/// running file works; on Windows the running `.exe` is locked, so it is moved
/// aside first (which Windows allows) and the new one renamed into its place.
fn replace_exe(new_bin: &Path, exe: &Path) -> std::io::Result<()> {
    if cfg!(windows) {
        let old = exe.with_extension("old");
        let _ = fs::remove_file(&old); // clear any leftover from a previous update
        fs::rename(exe, &old)?;
        fs::rename(new_bin, exe)?;
        let _ = fs::remove_file(&old); // best-effort; the running .old stays locked
        Ok(())
    } else {
        fs::rename(new_bin, exe)
    }
}

/// Delete the running executable. Unix unlinks it directly — the directory entry
/// goes away at once, the inode lives until this process exits. Windows can't
/// delete a locked running `.exe`, so it is moved aside (which Windows allows) and
/// a detached command deletes it once this process exits and releases the lock.
fn remove_exe(exe: &Path) -> std::io::Result<()> {
    if cfg!(windows) {
        let old = exe.with_extension("old");
        let _ = fs::remove_file(&old); // clear any leftover from a previous update
        fs::rename(exe, &old)?;
        schedule_self_delete(&old);
        Ok(())
    } else {
        fs::remove_file(exe)
    }
}

#[cfg(windows)]
fn schedule_self_delete(path: &Path) {
    use std::os::windows::process::CommandExt;
    // Detached, window-less cmd: pause for this process to exit (releasing the lock
    // on the moved-aside binary), then delete it. Best-effort — a leftover `.old`
    // is harmless and removable by hand.
    const CREATE_NO_WINDOW: u32 = 0x0800_0000;
    let line = format!(
        "/C ping 127.0.0.1 -n 2 >nul & del /F /Q \"{}\"",
        path.display()
    );
    let _ = Command::new("cmd")
        .raw_arg(line)
        .creation_flags(CREATE_NO_WINDOW)
        .spawn();
}

#[cfg(not(windows))]
fn schedule_self_delete(_path: &Path) {}

/// Removes its paths (files or directories) when dropped, so a failed update
/// leaves nothing behind.
struct CleanUp(Vec<PathBuf>);

impl Drop for CleanUp {
    fn drop(&mut self) {
        for p in &self.0 {
            let _ = fs::remove_file(p);
            let _ = fs::remove_dir_all(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newer_version_is_detected() {
        assert!(is_newer("0.2.0", "0.1.1"));
        assert!(is_newer("0.1.2", "0.1.1"));
        assert!(is_newer("1.0.0", "0.9.9"));
        assert!(!is_newer("0.1.1", "0.1.1"));
        assert!(!is_newer("0.1.0", "0.1.1"));
    }

    #[test]
    fn json_field_reads_the_tag() {
        let body = r#"{"url":"x","tag_name":"v0.2.3","name":"r"}"#;
        assert_eq!(json_field(body, "tag_name").as_deref(), Some("v0.2.3"));
        assert_eq!(json_field(body, "missing"), None);
    }

    #[test]
    fn semver_strips_v_and_suffix() {
        assert_eq!(semver("v1.2.3"), (1, 2, 3));
        assert_eq!(semver("1.2.3-rc1"), (1, 2, 3));
        assert_eq!(semver("0.1"), (0, 1, 0));
    }

    // Windows would move the file aside and spawn a deleter; exercise the direct
    // unlink that the `uninstall` subcommand uses on the user's platforms.
    #[cfg(unix)]
    #[test]
    fn remove_exe_deletes_the_binary() {
        let path =
            std::env::temp_dir().join(format!("isolf-uninstall-test-{}", std::process::id()));
        fs::write(&path, b"binary").unwrap();
        assert!(path.exists());
        remove_exe(&path).unwrap();
        assert!(!path.exists());
    }
}
