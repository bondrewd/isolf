//! Terminal UX: a small global reporter for sectioned, status-marked output.
//!
//! The build prints in phases (Build / Write), one status row each, a spinner on
//! the slow build step, and a one-line summary with the box, particle count,
//! elapsed time, and the next command. It adapts: colour and
//! spinners only on a terminal, ASCII glyphs on request, and `--quiet` /
//! `--verbose` levels. State lives in a process-global so the build code can
//! report from anywhere without threading a handle through every function.

use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, MutexGuard, OnceLock};
use std::time::{Duration, Instant};

use indicatif::{ProgressBar, ProgressStyle};

static STATE: OnceLock<Mutex<State>> = OnceLock::new();

struct State {
    color: bool,
    depth: Depth,
    unicode: bool,
    quiet: bool,
    verbose: bool,
    start: Instant,
    files: Vec<PathBuf>,
}

/// How much colour the terminal can render, so truecolour gradients degrade
/// gracefully instead of printing garbage.
#[derive(Clone, Copy, PartialEq)]
enum Depth {
    None,
    Ansi16,
    Ansi256,
    True,
}

/// Detect the colour depth from the environment (only when colour is wanted):
/// `COLORTERM=truecolor` → 24-bit, a `256` in `TERM` → 256-colour, else 16.
fn detect_depth(color: bool) -> Depth {
    if !color {
        return Depth::None;
    }
    let colorterm = std::env::var("COLORTERM").unwrap_or_default();
    if colorterm.contains("truecolor") || colorterm.contains("24bit") {
        return Depth::True;
    }
    if std::env::var("TERM").unwrap_or_default().contains("256") {
        return Depth::Ansi256;
    }
    Depth::Ansi16
}

/// Initialise the reporter. Call once at startup, before any other call.
pub fn init(color: bool, unicode: bool, quiet: bool, verbose: bool) {
    let _ = STATE.set(Mutex::new(State {
        color,
        depth: detect_depth(color),
        unicode,
        quiet,
        verbose,
        start: Instant::now(),
        files: Vec::new(),
    }));
}

fn state() -> MutexGuard<'static, State> {
    STATE
        .get()
        .expect("report::init must be called first")
        .lock()
        .unwrap()
}

// ---- colour and glyphs ----

/// Wrap `s` in the SGR `code` (e.g. "32" green, "1" bold, "2" dim) when colour is on.
fn paint(s: &str, code: &str, color: bool) -> String {
    if color && !code.is_empty() {
        format!("\x1b[{code}m{s}\x1b[0m")
    } else {
        s.to_string()
    }
}

/// Paint `s` in a 24-bit RGB colour, downsampled to what the terminal supports
/// (truecolour → 256 → 16 → no colour). `bold` adds the bold attribute.
fn rgb(s: &str, (r, g, b): (u8, u8, u8), depth: Depth, bold: bool) -> String {
    let pre = if bold { "1;" } else { "" };
    match depth {
        Depth::None => s.to_string(),
        Depth::True => format!("\x1b[{pre}38;2;{r};{g};{b}m{s}\x1b[0m"),
        Depth::Ansi256 => format!("\x1b[{pre}38;5;{}m{s}\x1b[0m", rgb_to_256(r, g, b)),
        Depth::Ansi16 => format!("\x1b[{pre}{}m{s}\x1b[0m", rgb_to_16(r, g, b)),
    }
}

/// Nearest xterm-256 index for an RGB colour (the 6×6×6 cube, or the grey ramp).
fn rgb_to_256(r: u8, g: u8, b: u8) -> u8 {
    let q = |v: u8| match v {
        0..=47 => 0,
        48..=114 => 1,
        115..=154 => 2,
        155..=194 => 3,
        195..=234 => 4,
        _ => 5,
    };
    if r == g && g == b {
        if r < 8 {
            16
        } else if r > 248 {
            231
        } else {
            232 + ((r as u16 - 8) * 24 / 247) as u8
        }
    } else {
        16 + 36 * q(r) + 6 * q(g) + q(b)
    }
}

/// Nearest of the 16 ANSI colours, as an SGR foreground code (30–37 / 90–97).
fn rgb_to_16(r: u8, g: u8, b: u8) -> u8 {
    let bright = r.max(g).max(b) > 160;
    let bit = |c: u8| u8::from(c > 110);
    let code = bit(r) + (bit(g) << 1) + (bit(b) << 2);
    (if bright { 90 } else { 30 }) + code
}

/// Linear interpolation between two RGB colours at `t` in `[0, 1]`.
fn lerp_rgb(a: (u8, u8, u8), b: (u8, u8, u8), t: f64) -> (u8, u8, u8) {
    let m = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * t).round() as u8;
    (m(a.0, b.0), m(a.1, b.1), m(a.2, b.2))
}

/// Display width of `s` in terminal cells, skipping ANSI SGR sequences so styled
/// rows measure by their visible text. The reporter's panel/title content is
/// single-width by construction (ASCII, box-drawing, block elements — all one
/// cell), so the visible character count is the cell count and borders stay square.
fn display_width(s: &str) -> usize {
    let mut width = 0;
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\x1b' {
            // Skip a CSI escape: ESC '[' … <final letter>.
            if chars.peek() == Some(&'[') {
                chars.next();
                for d in chars.by_ref() {
                    if d.is_ascii_alphabetic() {
                        break;
                    }
                }
            }
            continue;
        }
        width += 1;
    }
    width
}

struct Glyphs {
    ok: &'static str,
    warn: &'static str,
    fail: &'static str,
    dot: &'static str,
    arrow: &'static str,
}

fn glyphs(unicode: bool) -> Glyphs {
    if unicode {
        Glyphs {
            ok: "✓",
            warn: "⚠",
            fail: "✗",
            dot: "·",
            arrow: "→",
        }
    } else {
        Glyphs {
            ok: "ok",
            warn: "!",
            fail: "x",
            dot: "-",
            arrow: "->",
        }
    }
}

/// The middle-dot separator (`·`, or `-` in ASCII mode) used between detail parts.
pub fn dot() -> &'static str {
    glyphs(state().unicode).dot
}

/// Group a number into thousands (1234567 -> "1,234,567").
pub fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let bytes = digits.as_bytes();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

// ---- output primitives ----

/// The iSoLF wordmark, drawn with a colour gradient swept across its columns. The
/// `S`/`L`/`F` are figlet "ANSI Shadow"; the `i` and `o` are hand-set lowercase in
/// the same shadowed-block style (the force field is "iSoLF"). Single-width
/// box-drawing/block glyphs throughout, so it stays aligned.
const WORDMARK: [&str; 6] = [
    "   ███████╗         ██╗     ███████╗",
    "██╗██╔════╝ ██████╗ ██║     ██╔════╝",
    "╚═╝███████╗██╔═══██╗██║     █████╗  ",
    "██╗╚════██║██║   ██║██║     ██╔══╝  ",
    "██║███████║╚██████╔╝███████╗██║     ",
    "╚═╝╚══════╝ ╚═════╝ ╚══════╝╚═╝     ",
];

/// Gradient stops for the wordmark: teal → blue → violet.
const GRADIENT: [(u8, u8, u8); 3] = [(45, 212, 191), (59, 130, 246), (139, 92, 246)];

/// Tint each non-space cell of `text` by its column position, a smooth
/// left-to-right sweep through [`GRADIENT`]. Plain text when colour is off.
fn gradient_line(text: &str, width: usize, depth: Depth) -> String {
    if depth == Depth::None {
        return text.to_string();
    }
    let mut out = String::new();
    for (i, ch) in text.chars().enumerate() {
        if ch == ' ' {
            out.push(' ');
            continue;
        }
        let t = if width > 1 {
            i as f64 / (width - 1) as f64
        } else {
            0.0
        };
        let color = if t < 0.5 {
            lerp_rgb(GRADIENT[0], GRADIENT[1], t * 2.0)
        } else {
            lerp_rgb(GRADIENT[1], GRADIENT[2], (t - 0.5) * 2.0)
        };
        out.push_str(&rgb(&ch.to_string(), color, depth, true));
    }
    out
}

/// The header: a blank line, then the iSoLF wordmark with a colour gradient and
/// the version/tagline beside it (or a plain `isolf <version>` line in ASCII
/// mode). Omitted when quiet.
pub fn header() {
    let s = state();
    if s.quiet {
        return;
    }
    // A blank line so the wordmark is not flush against the shell's command line.
    println!();
    let version = concat!("v", env!("CARGO_PKG_VERSION"));
    if s.unicode {
        let width = WORDMARK.iter().map(|l| display_width(l)).max().unwrap_or(0);
        for (i, line) in WORDMARK.iter().enumerate() {
            let suffix = match i {
                2 => format!("   {}", paint(version, "2", s.color)),
                4 => format!("   {}", paint("coarse-grained", "2", s.color)),
                5 => format!("   {}", paint("membranes · vesicles", "2", s.color)),
                _ => String::new(),
            };
            println!("  {}{suffix}", gradient_line(line, width, s.depth));
        }
        println!();
    } else {
        println!(
            "{} {}\n",
            paint("isolf", "1", s.color),
            paint(version, "2", s.color)
        );
    }
}

/// A phase header (Build / Write …).
pub fn section(name: &str) {
    let s = state();
    if !s.quiet {
        println!("  {}", paint(name, "1", s.color));
    }
}

/// A success row: `✓ <label>  <detail>`. An empty label omits the column.
pub fn ok(label: &str, detail: &str) {
    let s = state();
    if s.quiet {
        return;
    }
    let mark = paint(glyphs(s.unicode).ok, "32", s.color);
    if label.is_empty() {
        println!("    {mark} {detail}");
    } else {
        let lab = paint(&format!("{label:<9}"), "1", s.color);
        println!("    {mark} {lab} {detail}");
    }
}

/// A warning row: `⚠ <msg>` (yellow). Shown even when quiet, since it matters.
pub fn warn(msg: &str) {
    let s = state();
    let mark = paint(glyphs(s.unicode).warn, "33", s.color);
    println!("    {mark} {}", paint(msg, "33", s.color));
}

/// Record a written file for the grouped Write section (see [`write_section`]).
pub fn record(path: &Path) {
    state().files.push(path.to_path_buf());
}

/// The files recorded so far (for the run log).
pub fn recorded_files() -> Vec<PathBuf> {
    state().files.clone()
}

/// Seconds since [`init`] (the build's elapsed time).
pub fn elapsed_secs() -> f64 {
    state().start.elapsed().as_secs_f64()
}

/// Render the Write section: every recorded file grouped by role under `dir`.
pub fn write_section(dir: &str) {
    let s = state();
    if s.quiet || s.files.is_empty() {
        return;
    }
    let g = glyphs(s.unicode);
    println!(
        "  {}  {}",
        paint("Write", "1", s.color),
        paint(dir, "36", s.color)
    );
    let roles: [(&str, &[&str]); 6] = [
        ("structure", &["gro", "pdb", "crd", "cif", "psf"]),
        ("topology", &["top", "itp"]),
        ("control", &["inp"]),
        ("view", &["vmd"]),
        ("animation", &["gif"]),
        ("log", &["log"]),
    ];
    for (role, exts) in roles {
        let names: Vec<String> = s
            .files
            .iter()
            .filter(|p| {
                p.extension()
                    .and_then(|e| e.to_str())
                    .is_some_and(|e| exts.contains(&e))
            })
            .map(|p| {
                if s.verbose {
                    p.display().to_string()
                } else {
                    p.file_name().unwrap_or_default().to_string_lossy().into()
                }
            })
            .collect();
        if !names.is_empty() {
            let label = paint(&format!("{role:<11}"), "2", s.color);
            println!("    {label} {}", names.join(&format!(" {} ", g.dot)));
        }
    }
}

/// Per-species colours for the composition bars (cycled by species index).
const SPECIES_COLORS: [(u8, u8, u8); 6] = [
    (96, 165, 250),  // blue
    (52, 211, 153),  // green
    (251, 191, 36),  // amber
    (244, 114, 182), // pink
    (167, 139, 250), // violet
    (45, 212, 191),  // teal
];

/// One leaflet for the composition panel: a label (`upper`/`lower` for a flat
/// membrane, `outer`/`inner` for a vesicle) and its species→count list, in the
/// membrane's species order.
pub type Leaflet = (&'static str, Vec<(String, usize)>);

/// The closing summary's inputs: identity, geometry, composition.
pub struct Summary<'a> {
    pub identity: &'a str,
    pub box_nm: [f64; 3],
    /// Lipid species and molecule counts (the composition panel).
    pub lipids: &'a [(String, usize)],
    /// Optional per-leaflet breakdown (the two leaflets, labelled), shown under the
    /// global composition.
    pub leaflets: Option<&'a [Leaflet]>,
    pub beads: usize,
    pub next: Option<String>,
}

/// Format the box dimensions, collapsing an equal-sided box to "L nm cube".
fn box_desc(box_nm: [f64; 3], unicode: bool) -> String {
    if (box_nm[0] - box_nm[1]).abs() < 1e-9 && (box_nm[1] - box_nm[2]).abs() < 1e-9 {
        format!("{:.1} nm cube", box_nm[0])
    } else {
        let x = if unicode { "×" } else { "x" };
        format!(
            "{:.0} {x} {:.0} {x} {:.0} nm",
            box_nm[0], box_nm[1], box_nm[2]
        )
    }
}

/// Teal accent for panel borders, and bright cyan for their titles.
const BORDER_ACCENT: (u8, u8, u8) = (45, 212, 191);
const TITLE_ACCENT: (u8, u8, u8) = (103, 232, 249);

/// Draw a rounded (or ASCII) panel: a bold bright-cyan title and teal borders,
/// with rows padded to a common width. Rows may carry colour; [`display_width`]
/// measures past it. The accent follows `depth`, so it downsamples and drops out.
fn panel(title: &str, rows: &[String], unicode: bool, depth: Depth) -> Vec<String> {
    let (tl, tr, bl, br, h, v) = if unicode {
        ("╭", "╮", "╰", "╯", "─", "│")
    } else {
        ("+", "+", "+", "+", "-", "|")
    };
    let inner = rows
        .iter()
        .map(|r| display_width(r))
        .max()
        .unwrap_or(0)
        .max(display_width(title) + 2);
    let span = inner + 2;
    // Top border: teal "corner dash space", cyan bold title, teal "space dashes corner".
    let trail = inner.saturating_sub(display_width(title) + 1);
    let mut out = vec![format!(
        "{}{}{}",
        rgb(&format!("{tl}{h} "), BORDER_ACCENT, depth, false),
        rgb(title, TITLE_ACCENT, depth, true),
        rgb(
            &format!(" {}{tr}", h.repeat(trail)),
            BORDER_ACCENT,
            depth,
            false
        ),
    )];
    let edge = rgb(v, BORDER_ACCENT, depth, false);
    for r in rows {
        let pad = inner - display_width(r);
        out.push(format!("{edge} {r}{} {edge}", " ".repeat(pad)));
    }
    out.push(rgb(
        &format!("{bl}{}{br}", h.repeat(span)),
        BORDER_ACCENT,
        depth,
        false,
    ));
    out
}

/// Place two panels side by side, padding the shorter to the taller's height.
fn hjoin(a: &[String], b: &[String], gap: usize) -> Vec<String> {
    let n = a.len().max(b.len());
    let wa = a.first().map_or(0, |l| display_width(l));
    let wb = b.first().map_or(0, |l| display_width(l));
    (0..n)
        .map(|i| {
            let la = a.get(i).cloned().unwrap_or_else(|| " ".repeat(wa));
            let lb = b.get(i).cloned().unwrap_or_else(|| " ".repeat(wb));
            format!("{la}{}{lb}", " ".repeat(gap))
        })
        .collect()
}

/// The closing summary: two side-by-side panels (composition, geometry) with the
/// next-command hint. Quiet mode collapses to a single status line instead.
pub fn summary(sum: &Summary) {
    let s = state();
    let g = glyphs(s.unicode);
    let total: usize = sum.lipids.iter().map(|(_, n)| n).sum();

    if s.quiet {
        let parts = [
            sum.identity.to_string(),
            box_desc(sum.box_nm, s.unicode),
            format!("{} beads", thousands(sum.beads)),
            format!("{:.1}s", s.start.elapsed().as_secs_f64()),
        ];
        println!(
            "  {} {}",
            paint(g.ok, "32", s.color),
            parts.join(&format!(" {} ", g.dot))
        );
        if let Some(next) = &sum.next {
            println!("    {} {next}", paint(g.arrow, "36", s.color));
        }
        return;
    }

    // Composition panel: one proportion bar per species, then a total.
    let (fill, empty) = if s.unicode {
        ("█", "░")
    } else {
        ("#", "-")
    };
    let mut comp = Vec::new();
    for (i, (name, count)) in sum.lipids.iter().enumerate() {
        let pct = (count * 100).checked_div(total).unwrap_or(0);
        let on = if total > 0 {
            (count * 10).div_ceil(total).min(10)
        } else {
            0
        };
        let bar = format!(
            "{}{}",
            rgb(&fill.repeat(on), SPECIES_COLORS[i % 6], s.depth, false),
            paint(&empty.repeat(10 - on), "2", s.color)
        );
        comp.push(format!("{name:<6} {count:>5}   {bar} {pct:>2}%"));
    }
    if !comp.is_empty() {
        comp.push(String::new());
    }
    comp.push(paint(
        &format!("total  {} lipids", thousands(total)),
        "2",
        s.color,
    ));
    // Per-leaflet breakdown: upper/lower for a membrane, outer/inner for a vesicle.
    if let Some(leaflets) = sum.leaflets {
        comp.push(String::new());
        // Right-align each count to the widest count in its column, so the dots
        // line up across the leaflet rows.
        let mut col_width: Vec<usize> = Vec::new();
        for (_, counts) in leaflets {
            for (j, (_, c)) in counts.iter().enumerate() {
                let w = thousands(*c).len();
                match col_width.get_mut(j) {
                    Some(cw) => *cw = (*cw).max(w),
                    None => col_width.push(w),
                }
            }
        }
        for (label, counts) in leaflets {
            let body = if counts.is_empty() {
                "—".to_string()
            } else {
                counts
                    .iter()
                    .enumerate()
                    .map(|(j, (sp, c))| format!("{sp} {:>w$}", thousands(*c), w = col_width[j]))
                    .collect::<Vec<_>>()
                    .join(" · ")
            };
            comp.push(format!(
                "{} {body}",
                paint(&format!("{label:<6}"), "2", s.color)
            ));
        }
    }
    let comp = panel("composition", &comp, s.unicode, s.depth);

    // Geometry panel: dim keys, plain values.
    let key = |k: &str| paint(&format!("{k:<8}"), "2", s.color);
    let mut geom = vec![
        format!("{} {}", key("box"), box_desc(sum.box_nm, s.unicode)),
        format!("{} {}", key("lipids"), thousands(total)),
    ];
    geom.push(format!("{} {}", key("beads"), thousands(sum.beads)));
    geom.push(format!(
        "{} {:.1} s",
        key("time"),
        s.start.elapsed().as_secs_f64()
    ));
    let geom = panel("geometry", &geom, s.unicode, s.depth);

    println!();
    for line in hjoin(&comp, &geom, 2) {
        println!("  {line}");
    }
    if let Some(next) = &sum.next {
        println!("  {} {next}", paint(g.arrow, "36", s.color));
    }
}

/// Print a fatal error (to stderr, red, with the ✗ marker). Safe before [`init`].
pub fn error(msg: &str) {
    let (color, unicode) = STATE
        .get()
        .map(|s| {
            let s = s.lock().unwrap();
            (s.color, s.unicode)
        })
        .unwrap_or((false, true));
    let mark = glyphs(unicode).fail;
    eprintln!("{}", paint(&format!("{mark} {msg}"), "1;31", color));
}

// ---- spinner ----

/// A live spinner over a slow step. Finishing it (explicitly or on drop) clears
/// the line; the caller then prints the matching `✓` row.
pub struct Spinner {
    bar: Option<ProgressBar>,
}

impl Spinner {
    /// Stop and clear the spinner now (otherwise it clears when dropped).
    pub fn finish(self) {}
}

impl Drop for Spinner {
    fn drop(&mut self) {
        if let Some(bar) = self.bar.take() {
            bar.finish_and_clear();
        }
    }
}

/// Start a spinner with `msg` on the slow step. No animation when quiet or when
/// stderr is not a terminal (the matching `✓` row still prints either way).
pub fn spin(msg: &str) -> Spinner {
    let s = state();
    if s.quiet || !std::io::stderr().is_terminal() {
        return Spinner { bar: None };
    }
    let ticks: &[&str] = if s.unicode {
        &["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
    } else {
        &["|", "/", "-", "\\"]
    };
    let bar = ProgressBar::new_spinner();
    bar.set_style(
        ProgressStyle::with_template("    {spinner:.cyan} {msg}")
            .unwrap()
            .tick_strings(ticks),
    );
    bar.set_message(msg.to_string());
    bar.enable_steady_tick(Duration::from_millis(80));
    Spinner { bar: Some(bar) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_width_skips_sgr_and_counts_cells() {
        assert_eq!(display_width("abc"), 3);
        // Colour codes do not add width; the box bars and borders are one cell each.
        assert_eq!(display_width("\x1b[1;38;2;1;2;3mabc\x1b[0m"), 3);
        assert_eq!(display_width("│ █████░░░░░ │"), 14);
    }

    #[test]
    fn panel_rows_are_all_equal_width() {
        // A panel pads every row (and its borders) to one width, so it stays square
        // even with a coloured row and a row wider than the title.
        let rows = [
            paint("short", "2", true),
            "a much longer row than the title".to_string(),
        ];
        // Truecolour borders/title carry SGR; display_width must still see one width.
        let p = panel("t", &rows, true, Depth::True);
        let widths: std::collections::BTreeSet<usize> =
            p.iter().map(|l| display_width(l)).collect();
        assert_eq!(widths.len(), 1, "panel rows differ in width: {p:?}");
    }

    #[test]
    fn hjoin_pads_the_shorter_panel() {
        let a = panel(
            "a",
            &["x".to_string(), "y".to_string(), "z".to_string()],
            true,
            Depth::None,
        );
        let b = panel("b", &["one".to_string()], true, Depth::None);
        let joined = hjoin(&a, &b, 2);
        assert_eq!(joined.len(), a.len()); // the taller panel's height
        let widths: std::collections::BTreeSet<usize> =
            joined.iter().map(|l| display_width(l)).collect();
        assert_eq!(widths.len(), 1, "joined rows differ in width");
    }

    #[test]
    fn rgb_downsampling_stays_in_range() {
        assert_eq!(rgb("x", (0, 0, 0), Depth::None, false), "x"); // no colour → untouched
        for &(r, g, b) in &[(0, 0, 0), (255, 255, 255), (45, 212, 191), (139, 92, 246)] {
            assert!((16..=255).contains(&rgb_to_256(r, g, b)));
            assert!((30..=97).contains(&rgb_to_16(r, g, b)));
        }
    }
}
