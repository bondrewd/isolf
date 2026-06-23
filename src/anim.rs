//! Animated GIF of the leaflet relaxation, written for `--gif`.
//!
//! Two render modes (`--gif-mode`):
//! - `point` (default): a top-down scatter, one disc per lipid. Clear for small to
//!   medium systems; a dense shell merges into a solid fill, so raise `--gif-scale`
//!   or switch to density there.
//! - `density`: a heatmap. Each lipid splats into an accumulation buffer and the
//!   pixels are coloured by local density, so a clumpy start (bright patches, dark
//!   gaps) smooths to an even mid-tone. Better for large systems.
//!
//! Both draw two panels (upper/outer left, lower/inner right). `--gif-scale`
//! multiplies the resolution and `--gif-fps` the playback speed. Frames are
//! subsampled to a cap, and the indexed buffers are LZW-encoded by the `gif` crate.

use std::fs::File;
use std::io;
use std::path::Path;

use isolf::membrane::RelaxFrame;

/// Pixels the longer box side maps to in each panel, before `--gif-scale`.
const PANEL: u32 = 360;
/// Border (px) around the panels.
const MARGIN: u32 = 16;
/// Gap (px) between the two panels.
const GAP: u32 = 24;
/// Most frames written; the relaxation is subsampled down to this.
const MAX_FRAMES: usize = 150;
/// Lipid disc radius in pixels, point mode (before `--gif-scale`).
const DOT_RADIUS: i32 = 2;
/// Per-lipid splat radius in pixels, density mode (before `--gif-scale`).
const SPLAT_RADIUS: i32 = 4;
/// Shades per leaflet in the density palette.
const RAMP: usize = 24;
/// Hold time for the final frame (~1.6 s) so the result is readable.
const HOLD_CS: u16 = 160;

/// Point-mode palette, RGB triples: background, upper leaflet, lower leaflet, field.
const PALETTE: &[u8] = &[
    0x1e, 0x1e, 0x2e, // 0 background
    0x2d, 0xd4, 0xbf, // 1 upper leaflet (teal)
    0x8b, 0x5c, 0xf6, // 2 lower leaflet (violet)
    0x31, 0x32, 0x44, // 3 panel field
];
const BG: u8 = 0;
const UPPER: u8 = 1;
const LOWER: u8 = 2;
const FIELD: u8 = 3;

/// How `--gif` draws the frames.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum GifMode {
    /// One disc per lipid (`--gif-mode point`).
    Point,
    /// Local-density heatmap (`--gif-mode density`).
    Density,
}

/// Settings from the `--gif-*` flags.
pub struct GifOptions {
    /// Render mode.
    pub mode: GifMode,
    /// Resolution multiplier (`--gif-scale`).
    pub scale: f64,
    /// Playback frames per second (`--gif-fps`).
    pub fps: u16,
}

/// Write `frames` as an animated GIF to `path`. `bounds` is the coordinate
/// rectangle `[min_x, min_y, max_x, max_y]` the recorded points live in (the box
/// for a membrane, the projection extent for a vesicle), used to scale them into
/// the panels.
pub fn write_gif(
    path: &Path,
    frames: &[RelaxFrame],
    bounds: [f64; 4],
    opts: &GifOptions,
) -> io::Result<()> {
    let scale = opts.scale.clamp(0.25, 8.0);
    let span_x = (bounds[2] - bounds[0]).max(1e-6);
    let span_y = (bounds[3] - bounds[1]).max(1e-6);
    let longer = span_x.max(span_y);
    let pw = ((span_x / longer) * PANEL as f64 * scale).round().max(1.0) as u32;
    let ph = ((span_y / longer) * PANEL as f64 * scale).round().max(1.0) as u32;
    let width = MARGIN + pw + GAP + pw + MARGIN;
    let height = MARGIN + ph + MARGIN;
    let splat_r = (SPLAT_RADIUS as f64 * scale).round().max(1.0) as i32;

    let palette: Vec<u8> = match opts.mode {
        GifMode::Point => PALETTE.to_vec(),
        GifMode::Density => density_palette(),
    };
    // GIF playback is a per-frame delay in centiseconds; turn fps into one.
    let delay = (100.0 / opts.fps.max(1) as f64).round().clamp(1.0, 255.0) as u16;
    let selected = subsample(frames.len(), MAX_FRAMES);

    let mut file = File::create(path)?;
    let mut encoder =
        gif::Encoder::new(&mut file, width as u16, height as u16, &palette).map_err(to_io)?;
    encoder.set_repeat(gif::Repeat::Infinite).map_err(to_io)?;

    for (n, &i) in selected.iter().enumerate() {
        let pixels = match opts.mode {
            GifMode::Point => draw_point(&frames[i], bounds, width, pw, ph),
            GifMode::Density => draw_density(&frames[i], bounds, width, pw, ph, splat_r),
        };
        let mut frame = gif::Frame::from_indexed_pixels(width as u16, height as u16, pixels, None);
        frame.delay = if n + 1 == selected.len() {
            HOLD_CS
        } else {
            delay
        };
        encoder.write_frame(&frame).map_err(to_io)?;
    }
    Ok(())
}

/// Indices of up to `cap` frames evenly spaced across `len`, always with the last.
fn subsample(len: usize, cap: usize) -> Vec<usize> {
    if len <= cap {
        return (0..len).collect();
    }
    let mut out: Vec<usize> = (0..cap).map(|k| k * (len - 1) / (cap - 1)).collect();
    out.dedup();
    out
}

/// Rasterise one frame as a scatter (point mode) into a `width * height` buffer.
fn draw_point(frame: &RelaxFrame, bounds: [f64; 4], width: u32, pw: u32, ph: u32) -> Vec<u8> {
    let height = MARGIN + ph + MARGIN;
    let span_x = (bounds[2] - bounds[0]).max(1e-6);
    let span_y = (bounds[3] - bounds[1]).max(1e-6);
    let mut buf = vec![BG; (width * height) as usize];
    let left = MARGIN;
    let right = MARGIN + pw + GAP;
    for (origin_x, color, points) in [(left, UPPER, &frame.upper), (right, LOWER, &frame.lower)] {
        fill_rect(&mut buf, width, origin_x, MARGIN, pw, ph, FIELD);
        for &[x, y] in points {
            // Map (x, y) within the bounds to panel pixels; flip y so +y points up.
            let px = origin_x as f64 + ((x - bounds[0]) / span_x).clamp(0.0, 1.0) * (pw - 1) as f64;
            let py = MARGIN as f64
                + (1.0 - ((y - bounds[1]) / span_y).clamp(0.0, 1.0)) * (ph - 1) as f64;
            draw_dot(&mut buf, width, height, px as i32, py as i32, color);
        }
    }
    buf
}

/// Rasterise one frame as a local-density heatmap (density mode).
fn draw_density(
    frame: &RelaxFrame,
    bounds: [f64; 4],
    width: u32,
    pw: u32,
    ph: u32,
    splat_r: i32,
) -> Vec<u8> {
    let height = MARGIN + ph + MARGIN;
    let span_x = (bounds[2] - bounds[0]).max(1e-6);
    let span_y = (bounds[3] - bounds[1]).max(1e-6);
    let mut buf = vec![BG; (width * height) as usize];
    let left = MARGIN;
    let right = MARGIN + pw + GAP;
    let outer_base = 1u8;
    let inner_base = 1 + RAMP as u8;
    for (origin_x, base, points) in [
        (left, outer_base, &frame.upper),
        (right, inner_base, &frame.lower),
    ] {
        // Splat every lipid into an accumulation buffer; overlapping splats build
        // up the local density.
        let mut accum = vec![0.0f32; (pw * ph) as usize];
        for &[x, y] in points {
            let px = (((x - bounds[0]) / span_x).clamp(0.0, 1.0) * (pw - 1) as f64).round() as i32;
            let py = ((1.0 - ((y - bounds[1]) / span_y).clamp(0.0, 1.0)) * (ph - 1) as f64).round()
                as i32;
            splat(&mut accum, pw, ph, px, py, splat_r);
        }
        // Calibrate against the mean over the covered area, so an even shell sits
        // mid-ramp: brighter where denser, darker where sparser.
        let sum: f64 = accum.iter().map(|&v| v as f64).sum();
        let covered = accum.iter().filter(|&&v| v > 0.0).count();
        let reference = if covered > 0 {
            sum / covered as f64 * 1.6
        } else {
            1.0
        };
        for y in 0..ph {
            for x in 0..pw {
                let v = accum[(y * pw + x) as usize] as f64 / reference;
                let idx = (v * (RAMP - 1) as f64)
                    .round()
                    .clamp(0.0, (RAMP - 1) as f64) as u8;
                buf[((MARGIN + y) * width + origin_x + x) as usize] = base + idx;
            }
        }
    }
    buf
}

/// Add a filled disc of weight to `accum`, centred at `(cx, cy)` and clipped to the panel.
fn splat(accum: &mut [f32], pw: u32, ph: u32, cx: i32, cy: i32, r: i32) {
    for dy in -r..=r {
        for dx in -r..=r {
            if dx * dx + dy * dy > r * r {
                continue;
            }
            let (x, y) = (cx + dx, cy + dy);
            if x >= 0 && y >= 0 && (x as u32) < pw && (y as u32) < ph {
                accum[(y as u32 * pw + x as u32) as usize] += 1.0;
            }
        }
    }
}

/// Density palette: background, then the outer (teal) and inner (violet) ramps,
/// each `RAMP` shades dark-to-bright, padded to a power-of-two colour count.
fn density_palette() -> Vec<u8> {
    let mut palette = vec![0x1e, 0x1e, 0x2e]; // 0 background
    for (low, high) in [
        ((0x24, 0x26, 0x38), (0x5e, 0xf0, 0xd8)), // outer: dark -> teal
        ((0x26, 0x22, 0x3a), (0xa9, 0x86, 0xff)), // inner: dark -> violet
    ] {
        for i in 0..RAMP {
            let (r, g, b) = lerp(low, high, i as f64 / (RAMP - 1) as f64);
            palette.extend_from_slice(&[r, g, b]);
        }
    }
    let colors = (palette.len() / 3).max(1);
    palette.resize(colors.next_power_of_two() * 3, 0);
    palette
}

/// Linearly interpolate between two RGB colours, `t` in `[0, 1]`.
fn lerp(a: (u8, u8, u8), b: (u8, u8, u8), t: f64) -> (u8, u8, u8) {
    let mix = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * t).round() as u8;
    (mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

/// Fill the `w` by `h` rectangle at `(x0, y0)` with `color`.
fn fill_rect(buf: &mut [u8], width: u32, x0: u32, y0: u32, w: u32, h: u32, color: u8) {
    for y in y0..y0 + h {
        let row = (y * width) as usize;
        buf[row + x0 as usize..row + (x0 + w) as usize].fill(color);
    }
}

/// Draw a filled disc of `color` centred at `(cx, cy)`, clipped to the canvas.
fn draw_dot(buf: &mut [u8], width: u32, height: u32, cx: i32, cy: i32, color: u8) {
    let r = DOT_RADIUS;
    for dy in -r..=r {
        for dx in -r..=r {
            if dx * dx + dy * dy > r * r {
                continue;
            }
            let (x, y) = (cx + dx, cy + dy);
            if x >= 0 && y >= 0 && (x as u32) < width && (y as u32) < height {
                buf[(y as u32 * width + x as u32) as usize] = color;
            }
        }
    }
}

/// Wrap a `gif` encoding error as an [`io::Error`].
fn to_io<E: std::fmt::Display>(error: E) -> io::Error {
    io::Error::other(error.to_string())
}
