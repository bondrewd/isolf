//! The always-on run log: a plain-text provenance file written next to the
//! output, recording the exact (re-runnable) command, the seed, who/where/when it
//! ran, and what it produced. Reading it tells you how to reproduce the run.
//!
//! No extra dependencies: the timestamp is formatted from `SystemTime` by hand
//! (UTC), and the user and host come from the environment with a `hostname`
//! fallback.

use std::fmt::Write as _;
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// One run's worth of provenance, ready to render to a log file.
pub struct Run<'a> {
    /// The exact, re-runnable command (seed pinned in).
    pub command: &'a str,
    /// The seed actually used.
    pub seed: u64,
    /// Build wall-clock time (seconds).
    pub duration_secs: f64,
    /// Build-mode label (e.g. "vesicle").
    pub identity: &'a str,
    /// Final system box (nm).
    pub box_nm: [f64; 3],
    /// Total particle count.
    pub beads: usize,
    /// Output directory, as given.
    pub out_dir: &'a str,
    /// File names written (the log itself excluded).
    pub files: &'a [String],
}

impl Run<'_> {
    /// Render and write the log to `path`.
    pub fn write(&self, path: &Path) -> std::io::Result<()> {
        std::fs::write(path, self.render())
    }

    fn render(&self) -> String {
        let mut s = String::new();
        let _ = writeln!(s, "isolf {} run log", env!("CARGO_PKG_VERSION"));
        let _ = writeln!(s);
        let _ = writeln!(s, "command    {}", self.command);
        let _ = writeln!(s, "date       {}", utc_now());
        let _ = writeln!(s, "host       {}", host());
        let _ = writeln!(s, "user       {}", user());
        let _ = writeln!(s, "cwd        {}", cwd());
        let _ = writeln!(s, "seed       {}", self.seed);
        let _ = writeln!(s, "duration   {:.1} s", self.duration_secs);
        let _ = writeln!(s);
        let _ = writeln!(s, "mode       {}", self.identity);
        let _ = writeln!(s, "box        {} nm", box_dims(self.box_nm));
        let _ = writeln!(s, "beads      {}", self.beads);
        let _ = writeln!(s, "output     {}", self.out_dir);
        let _ = writeln!(s, "files      {}", self.files.join(", "));
        s
    }
}

/// "40.00" if cubic, else "40.00 x 40.00 x 15.62" (ASCII, so the log greps clean).
fn box_dims(b: [f64; 3]) -> String {
    if (b[0] - b[1]).abs() < 1e-9 && (b[1] - b[2]).abs() < 1e-9 {
        format!("{:.2}", b[0])
    } else {
        format!("{:.2} x {:.2} x {:.2}", b[0], b[1], b[2])
    }
}

/// The current UTC time as `YYYY-MM-DD HH:MM:SS UTC`, from the Unix epoch with
/// Howard Hinnant's civil-from-days conversion (no calendar dependency).
fn utc_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    let days = secs.div_euclid(86_400);
    let tod = secs.rem_euclid(86_400);
    let (hour, min, sec) = (tod / 3600, (tod % 3600) / 60, tod % 60);

    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097; // day of era [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365; // [0, 399]
    let year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // day of year [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let day = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let month = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    let year = if month <= 2 { year + 1 } else { year };

    format!("{year:04}-{month:02}-{day:02} {hour:02}:{min:02}:{sec:02} UTC")
}

/// The hostname: `$HOSTNAME` / `$COMPUTERNAME`, else the `hostname` command, else
/// "unknown".
fn host() -> String {
    if let Some(h) = env_first(&["HOSTNAME", "COMPUTERNAME"]) {
        return h;
    }
    if let Ok(out) = Command::new("hostname").output()
        && out.status.success()
    {
        let name = String::from_utf8_lossy(&out.stdout).trim().to_string();
        if !name.is_empty() {
            return name;
        }
    }
    "unknown".to_string()
}

/// The user: `$USER` / `$LOGNAME` / `$USERNAME`, else "unknown".
fn user() -> String {
    env_first(&["USER", "LOGNAME", "USERNAME"]).unwrap_or_else(|| "unknown".to_string())
}

fn cwd() -> String {
    std::env::current_dir()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| "unknown".to_string())
}

fn env_first(keys: &[&str]) -> Option<String> {
    keys.iter()
        .filter_map(|k| std::env::var(k).ok())
        .find(|v| !v.trim().is_empty())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn utc_known_epoch() {
        // 1_700_000_000 = 2023-11-14 22:13:20 UTC (a fixed reference).
        let secs = 1_700_000_000_i64;
        let days = secs.div_euclid(86_400);
        let z = days + 719_468;
        let era = z.div_euclid(146_097);
        let doe = z - era * 146_097;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let year = yoe + era * 400;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let day = doy - (153 * mp + 2) / 5 + 1;
        let month = if mp < 10 { mp + 3 } else { mp - 9 };
        let year = if month <= 2 { year + 1 } else { year };
        assert_eq!((year, month, day), (2023, 11, 14));
    }

    #[test]
    fn box_dims_collapses_cubic() {
        assert_eq!(box_dims([40.0, 40.0, 40.0]), "40.00");
        assert_eq!(box_dims([40.0, 40.0, 15.62]), "40.00 x 40.00 x 15.62");
    }

    #[test]
    fn render_has_the_essentials() {
        let run = Run {
            command: "isolf --upper POPC=1 --membrane 40 --seed 7 --out m",
            seed: 7,
            duration_secs: 2.4,
            identity: "membrane",
            box_nm: [40.0, 40.0, 15.62],
            beads: 6348,
            out_dir: "m",
            files: &["m.gro".to_string(), "m.top".to_string()],
        };
        let text = run.render();
        assert!(text.contains("seed       7"));
        assert!(text.contains("--seed 7"));
        assert!(text.contains("command    isolf"));
        assert!(text.contains("m.gro, m.top"));
    }
}
