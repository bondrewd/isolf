//! Construction of a coarse-grained lipid bilayer.
//!
//! Each lipid is dropped at a uniform-random `(x, y)` in the shared periodic box
//! and laid out along z by walking its bond lengths outward from the phosphate
//! bead, so head groups point away from the bilayer midplane and tails toward
//! it. The upper leaflet keeps this orientation; the lower leaflet is mirrored
//! through z. The random start is then relaxed to an even, clash-free layout: a
//! per-leaflet Lloyd (centroidal-Voronoi) pass and then a soft-sphere de-clash,
//! both in [`crate::pack`].
//!
//! **Lipids start at the clash-free minimal-energy spacing** (`√1.08 · 2^(1/6) ·
//! σ_max`), just past the bead potential minimum: a low-energy start that an NPT
//! (semi-isotropic, Z-fixed) run then compresses to the real density. The box
//! holds each leaflet's lipid count at that spacing, sized as a square to
//! whichever leaflet needs more area in count mode, or to the x/y sides given in
//! box mode. The lighter leaflet sits at lower density until NPT compresses the
//! shared box.

use std::f64::consts::TAU;

use rand::Rng;
use rayon::prelude::*;

use crate::composition::Composition;
use crate::error::BuildError;
use crate::force_field::{ForceField, POTENTIAL_MINIMUM_FACTOR};
use crate::pack;
use crate::voronoi::wrap;

/// Default system name written into the output files.
pub const DEFAULT_NAME: &str = "CG membrane model";

/// Default system temperature (K): the default temperature of the generated
/// GENESIS control files.
pub const DEFAULT_TEMPERATURE: f64 = 323.15;

/// Default vacuum padding (nm) added on each side: along z for a flat membrane,
/// on every side for a vesicle.
pub const DEFAULT_PADDING: f64 = 10.0;

/// Area expansion over close packing for the default energy-minimizing layout:
/// lipids sit just past the bead potential minimum, so the start is low-energy
/// and clash-free (the small extra area is what an NPT run then compresses out).
pub(crate) const ENERGY_PACKING_AREA_FACTOR: f64 = 1.08;

/// Default relaxation budget for packing lipids: iterate until the largest bead
/// overlap drops below the tolerance (nm) or the iteration cap is hit.
pub const DEFAULT_PACK_ITERATIONS: usize = 1000;
pub const DEFAULT_PACK_TOLERANCE: f64 = 0.01;

/// A single coarse-grained particle (one bead of one lipid) placed in space.
#[derive(Debug, Clone, PartialEq)]
pub struct Particle {
    /// Lipid residue name (e.g. `POPC`).
    pub residue: String,
    /// Bead-type name (e.g. `PHO`); also serves as the atom type.
    pub bead: String,
    /// 1-based molecule (lipid) index, continuous across both leaflets.
    pub molecule_id: usize,
    /// 1-based atom (bead) index, continuous across both leaflets.
    pub atom_id: usize,
    /// Cartesian position (nm).
    pub position: [f64; 3],
    /// Bead charge (e), copied from the force field so each output format can be
    /// written from the [`Membrane`] alone.
    pub charge: f64,
    /// Bead mass (g/mol), copied from the force field.
    pub mass: f64,
}

/// How to size the membrane in-plane: by a lipid count per leaflet (a square
/// box), by an explicit box (`x`, `y` in nm; rectangular when `x != y`), or by a
/// count placed into a given box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Sizing {
    /// Pack the given number of lipids per leaflet; the builder picks the smallest
    /// square box that holds them at the clash-free minimal-energy spacing.
    Count { upper: usize, lower: usize },
    /// Fill a box of the given in-plane size (nm) with as many lipids as fit per
    /// leaflet. `x == y` is a square membrane.
    Box { x: f64, y: f64 },
    /// Pack the given number of lipids per leaflet into a box of the given
    /// in-plane size (nm).
    CountInBox {
        upper: usize,
        lower: usize,
        x: f64,
        y: f64,
    },
}

/// Lateral lipid organization within each leaflet or shell.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Phase {
    /// Uniform random mixing across the entire leaflet/shell (default).
    #[default]
    Random,
    /// Contiguous stripes (flat membrane), latitude bands (vesicle), or axial
    /// rings (tubule), sized to each domain's area fraction.
    Stripe,
    /// Circular raft (flat membrane / tubule) or polar cap (vesicle) formed by
    /// the first domain, surrounded by the remaining lipids.
    Raft,
}

/// Knobs that tune the build independently of the sizing mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BuildOptions {
    /// Vacuum padding (nm) added on each side: along z for a flat membrane, on
    /// every side for a vesicle, and along x/y for a periodic tubule.
    pub padding: f64,
    /// Lateral domain layout within each leaflet/shell.
    pub phase: Phase,
    /// Offset phase domains between the two leaflets/shells (half-box shift on
    /// membranes, opposite pole on vesicles, half-length + π shift on tubules).
    pub anti_register: bool,
}

impl Default for BuildOptions {
    fn default() -> Self {
        Self {
            padding: DEFAULT_PADDING,
            phase: Phase::Random,
            anti_register: false,
        }
    }
}

/// A built membrane: every particle plus the metadata the writers need.
#[derive(Debug, Clone, PartialEq)]
pub struct Membrane {
    /// System name (used in output file titles).
    pub name: String,
    /// All particles, in the canonical layout shared by every writer: grouped by
    /// lipid species (matching [`lipid_counts`](Membrane::lipid_counts)) with
    /// molecule and atom ids numbered sequentially from one.
    pub particles: Vec<Particle>,
    /// Periodic box (nm): the in-plane box (square, or the given x/y sides), with
    /// z the membrane thickness plus padding on each side.
    pub box_size: [f64; 3],
    /// Molecule count per lipid species, in a stable canonical order shared by
    /// every output file.
    pub lipid_counts: Vec<(String, usize)>,
}

/// One leaflet's lipid `count` and the `cells_x · cells_y` cell lattice (at the
/// lipid spacing) that sizes its share of the box. Lipids are placed at random
/// within the box, not on this lattice.
#[derive(Debug, Clone, Copy)]
struct LeafletGrid {
    cells_x: usize,
    cells_y: usize,
    count: usize,
}

impl LeafletGrid {
    fn capacity(self) -> usize {
        self.cells_x * self.cells_y
    }
}

/// The resolved geometry: a shared in-plane box (nm) and a grid per leaflet.
#[derive(Debug, Clone, Copy)]
struct Plan {
    box_x: f64,
    box_y: f64,
    upper: LeafletGrid,
    lower: LeafletGrid,
}

/// One captured step of the lateral relaxation, for the `--gif` animation: the
/// `(x, y)` position of each lipid in each leaflet (one point per lipid), in the
/// build's pre-centred box `[0, box]`.
#[derive(Debug, Clone)]
pub struct RelaxFrame {
    /// Upper-leaflet lipid positions (nm).
    pub upper: Vec<[f64; 2]>,
    /// Lower-leaflet lipid positions (nm).
    pub lower: Vec<[f64; 2]>,
}

impl Membrane {
    /// Build a bilayer from upper and lower leaflet compositions.
    ///
    /// `rng` seeds each lipid's random position; the relaxation that follows is
    /// deterministic, so a seeded generator reproduces the output.
    ///
    /// # Errors
    ///
    /// Returns a [`BuildError`] if a composition names an unknown lipid or bead, a
    /// lipid has no phosphate bead, a leaflet ends up empty, or more lipids are
    /// requested than the box holds.
    pub fn build(
        force_field: &ForceField,
        name: &str,
        upper: &Composition,
        lower: &Composition,
        sizing: Sizing,
        options: &BuildOptions,
        rng: &mut impl Rng,
    ) -> Result<Self, BuildError> {
        Self::build_inner(force_field, name, upper, lower, sizing, options, rng, None)
    }

    /// Like [`build`](Membrane::build), but also records every relaxation step (the
    /// Lloyd rounds, then the de-clash) as a [`RelaxFrame`], for the `--gif`
    /// animation.
    ///
    /// # Errors
    ///
    /// The same as [`build`](Membrane::build).
    pub fn build_recorded(
        force_field: &ForceField,
        name: &str,
        upper: &Composition,
        lower: &Composition,
        sizing: Sizing,
        options: &BuildOptions,
        rng: &mut impl Rng,
    ) -> Result<(Self, Vec<RelaxFrame>), BuildError> {
        let mut frames = Vec::new();
        let membrane = Self::build_inner(
            force_field,
            name,
            upper,
            lower,
            sizing,
            options,
            rng,
            Some(&mut frames),
        )?;
        Ok((membrane, frames))
    }

    #[allow(clippy::too_many_arguments)]
    fn build_inner(
        force_field: &ForceField,
        name: &str,
        upper: &Composition,
        lower: &Composition,
        sizing: Sizing,
        options: &BuildOptions,
        rng: &mut impl Rng,
        frames: Option<&mut Vec<RelaxFrame>>,
    ) -> Result<Self, BuildError> {
        let upper_spacing = leaflet_spacing(force_field, upper)?;
        let lower_spacing = leaflet_spacing(force_field, lower)?;
        let plan = resolve_geometry(sizing, upper_spacing, lower_spacing)?;

        let mut molecule_id = 0;
        let mut atom_id = 0;
        let mut particles = Vec::new();
        let mut radii = Vec::new();
        let mut ranges = Vec::new();

        let upper_counts = place_leaflet(
            force_field,
            upper,
            plan.box_x,
            plan.box_y,
            plan.upper,
            upper_spacing,
            Leaflet::Upper,
            options,
            rng,
            &mut molecule_id,
            &mut atom_id,
            &mut particles,
            &mut radii,
            &mut ranges,
        )?;
        // The lipids placed so far form the upper leaflet (the lower follows).
        let upper_lipids = molecule_id;
        let lower_counts = place_leaflet(
            force_field,
            lower,
            plan.box_x,
            plan.box_y,
            plan.lower,
            lower_spacing,
            Leaflet::Lower,
            options,
            rng,
            &mut molecule_id,
            &mut atom_id,
            &mut particles,
            &mut radii,
            &mut ranges,
        )?;

        // Relax the random placement to a uniform, clash-free start: an independent
        // centroidal-Voronoi (Lloyd) pass per leaflet, then a shared de-clash.
        relax_membrane(
            &mut particles,
            &radii,
            &ranges,
            upper_lipids,
            plan.box_x,
            plan.box_y,
            frames,
        );

        // Centre the bilayer on the origin in x/y (z is already on the midplane), the
        // convention the output writers expect.
        for particle in &mut particles {
            particle.position[0] -= plan.box_x / 2.0;
            particle.position[1] -= plan.box_y / 2.0;
        }

        let lipid_counts = aggregate_counts(upper, &upper_counts, lower, &lower_counts);
        // Reorder into one canonical layout (grouped by species) shared by every
        // output file, so any coordinate file matches any topology file.
        let particles = canonicalize(particles, &lipid_counts);

        Ok(Self {
            name: name.to_string(),
            box_size: periodic_box(&particles, plan.box_x, plan.box_y, options.padding),
            lipid_counts,
            particles,
        })
    }

    /// Total number of lipid molecules in the membrane.
    pub fn total_lipids(&self) -> usize {
        self.lipid_counts.iter().map(|(_, count)| count).sum()
    }

    /// Translate every particle in place so the particle centroid sits at the
    /// centre of the box `[0, box_size]` on each axis.
    pub fn center_in_box(&mut self) {
        let n = self.particles.len() as f64;
        if n == 0.0 {
            return;
        }
        let mut sum = [0.0; 3];
        for p in &self.particles {
            sum[0] += p.position[0];
            sum[1] += p.position[1];
            sum[2] += p.position[2];
        }
        self.translate([
            self.box_size[0] / 2.0 - sum[0] / n,
            self.box_size[1] / 2.0 - sum[1] / n,
            self.box_size[2] / 2.0 - sum[2] / n,
        ]);
    }

    /// A copy translated so the particle centroid sits at the centre of the box
    /// `[0, box_size]` on each axis. The build leaves a flat bilayer centred on the
    /// origin (so it draws in a corner of the periodic box a viewer renders); this
    /// places it centred, matching the vesicle path. It is a no-op for an
    /// already-centred membrane, e.g. a vesicle.
    pub fn centered_in_box(&self) -> Membrane {
        let mut centered = self.clone();
        centered.center_in_box();
        centered
    }

    /// Translate every particle by `shift` (nm), in place. Used to re-centre a
    /// finished build on the coordinate origin for `--center origin`.
    pub fn translate(&mut self, shift: [f64; 3]) {
        for p in &mut self.particles {
            for (pos, s) in p.position.iter_mut().zip(shift) {
                *pos += s;
            }
        }
    }
}

/// Which leaflet is being placed: the lower one is mirrored through z.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Leaflet {
    Upper,
    Lower,
}

/// Lateral spacing (nm) for a leaflet: the clash-free minimal-energy spacing
/// `√1.08 · 2^(1/6)·σ_max`, placing lipids just past the bead potential minimum
/// (a low-energy start that an NPT run compresses to the real density).
pub(crate) fn leaflet_spacing(
    force_field: &ForceField,
    composition: &Composition,
) -> Result<f64, BuildError> {
    Ok(ENERGY_PACKING_AREA_FACTOR.sqrt() * clash_floor(force_field, composition)?)
}

/// The lateral area one lipid occupies in a leaflet (nm²): the square of the
/// lipid spacing. Used to match the two leaflets' areas.
///
/// # Errors
///
/// Returns a [`BuildError`] if the composition is empty or names a lipid or bead
/// the force field does not define.
pub fn area_per_lipid(
    force_field: &ForceField,
    composition: &Composition,
) -> Result<f64, BuildError> {
    let spacing = leaflet_spacing(force_field, composition)?;
    Ok(spacing * spacing)
}

/// The clash-free lower bound on the spacing: `2^(1/6) · σ_max` over the beads of
/// the lipids present (the densest a leaflet may start without bead overlap).
pub(crate) fn clash_floor(
    force_field: &ForceField,
    composition: &Composition,
) -> Result<f64, BuildError> {
    let mut minimum = 0.0_f64;
    for component in composition.components() {
        let lipid = force_field
            .lipid(&component.name)
            .ok_or_else(|| BuildError::UnknownLipid(component.name.clone()))?;
        for bead_name in &lipid.beads {
            let bead = force_field
                .bead_type(bead_name)
                .ok_or_else(|| BuildError::UnknownBead(bead_name.clone()))?;
            minimum = minimum.max(POTENTIAL_MINIMUM_FACTOR * bead.interaction.sigma());
        }
    }
    if minimum <= 0.0 {
        return Err(BuildError::EmptyComposition);
    }
    Ok(minimum)
}

/// Resolve the shared in-plane box and per-leaflet grids for the chosen sizing,
/// given each leaflet's lateral spacing.
fn resolve_geometry(
    sizing: Sizing,
    upper_spacing: f64,
    lower_spacing: f64,
) -> Result<Plan, BuildError> {
    match sizing {
        // Smallest square box that holds each leaflet's count at its spacing.
        Sizing::Count { upper, lower } => {
            let up = square_grid_for_count(upper);
            let lo = square_grid_for_count(lower);
            let side = (up.cells_x as f64 * upper_spacing).max(lo.cells_x as f64 * lower_spacing);
            if side <= 0.0 {
                return Err(BuildError::EmptyLeaflet);
            }
            Ok(Plan {
                box_x: side,
                box_y: side,
                upper: up,
                lower: lo,
            })
        }
        // Fill the given box with as many lipids as fit at each leaflet's spacing.
        Sizing::Box { x, y } => {
            let up = fill_grid_for_box(x, y, upper_spacing)?;
            let lo = fill_grid_for_box(x, y, lower_spacing)?;
            Ok(Plan {
                box_x: x,
                box_y: y,
                upper: up,
                lower: lo,
            })
        }
        // Pack the requested counts into the given box.
        Sizing::CountInBox { upper, lower, x, y } => {
            let up = pack_grid_into_box(x, y, upper_spacing, upper)?;
            let lo = pack_grid_into_box(x, y, lower_spacing, lower)?;
            Ok(Plan {
                box_x: x,
                box_y: y,
                upper: up,
                lower: lo,
            })
        }
    }
}

/// Smallest square grid holding `count` lipids (empty when `count` is 0).
fn square_grid_for_count(count: usize) -> LeafletGrid {
    let n = (count as f64).sqrt().ceil() as usize;
    LeafletGrid {
        cells_x: n,
        cells_y: n,
        count,
    }
}

/// Grid that fills an `x` by `y` box at `spacing`, one lipid per cell.
fn fill_grid_for_box(x: f64, y: f64, spacing: f64) -> Result<LeafletGrid, BuildError> {
    let cells_x = (x / spacing).floor() as usize;
    let cells_y = (y / spacing).floor() as usize;
    if cells_x == 0 || cells_y == 0 {
        return Err(BuildError::EmptyLeaflet);
    }
    Ok(LeafletGrid {
        cells_x,
        cells_y,
        count: cells_x * cells_y,
    })
}

/// Grid that packs `count` lipids into an `x` by `y` box at `spacing`.
fn pack_grid_into_box(
    x: f64,
    y: f64,
    spacing: f64,
    count: usize,
) -> Result<LeafletGrid, BuildError> {
    let cells_x = (x / spacing).floor() as usize;
    let cells_y = (y / spacing).floor() as usize;
    let grid = LeafletGrid {
        cells_x,
        cells_y,
        count,
    };
    if count > grid.capacity() {
        return Err(BuildError::BoxTooSmall {
            requested: count,
            capacity: grid.capacity(),
        });
    }
    Ok(grid)
}

/// Area share of each phase domain (`0..composition.domain_count()`) weighted
/// by each species' packing area (`(2^(1/6) · σ_max)²`), so every domain starts
/// at the same local areal density and shared Lloyd CVT preserves domain walls.
pub(crate) fn domain_area_fractions(
    composition: &Composition,
    counts: &[usize],
    force_field: &ForceField,
) -> Vec<f64> {
    let mut areas = vec![0.0_f64; composition.domain_count()];
    for (idx, (component, &count)) in composition.components().iter().zip(counts).enumerate() {
        if count == 0 {
            continue;
        }
        let max_sigma = force_field
            .lipid(&component.name)
            .map(|l| {
                l.beads
                    .iter()
                    .filter_map(|b| force_field.bead_type(b))
                    .map(|bt| bt.interaction.sigma())
                    .fold(0.0_f64, f64::max)
            })
            .unwrap_or(1.0);
        let d = POTENTIAL_MINIMUM_FACTOR * max_sigma;
        areas[composition.domain_of(idx)] += count as f64 * d * d;
    }
    let total: f64 = areas.iter().sum();
    if total > 0.0 {
        for area in &mut areas {
            *area /= total;
        }
    }
    areas
}

/// Exact area of the disk `{x² + y² ≤ r²}` clipped to `[-a, a] × [-b, b]`.
fn disk_in_rect_area(r: f64, a: f64, b: f64) -> f64 {
    if r <= 0.0 || a <= 0.0 || b <= 0.0 {
        return 0.0;
    }
    let r2 = r * r;
    if r2 >= a * a + b * b {
        return 4.0 * a * b;
    }
    let x_hi = a.min(r);
    let x_flat = if b < r {
        (r2 - b * b).sqrt().min(x_hi)
    } else {
        0.0
    };
    let antideriv = |x: f64| {
        let y = (r2 - x * x).max(0.0).sqrt();
        0.5 * (x * y + r2 * (x / r).clamp(-1.0, 1.0).asin())
    };
    let quarter = b * x_flat + (antideriv(x_hi) - antideriv(x_flat));
    4.0 * quarter
}

/// Solve for the radius `r` such that a disk clipped to `[0, bx) × [0, by)`
/// covers exact area fraction `fraction ∈ [0, 1]`.
fn solve_raft_radius(bx: f64, by: f64, fraction: f64) -> f64 {
    let a = 0.5 * bx;
    let b = 0.5 * by;
    let target = fraction.clamp(0.0, 1.0) * bx * by;
    let mut lo = 0.0_f64;
    let mut hi = a.hypot(b);
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        if disk_in_rect_area(mid, a, b) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

enum PlaneDomainMode {
    Random,
    Stripe {
        axis: usize,
        shift: f64,
        bounds: Vec<(f64, f64)>,
    },
    Raft {
        center: [f64; 2],
        radius: f64,
    },
}

/// Phase-domain sampler on a periodic 2D rectangle `[0, box_xy[0]) × [0, box_xy[1])`.
/// Shared by flat bilayers (`stripe_axis = 0`) and unrolled cylindrical tubules
/// (`stripe_axis = 1`).
pub(crate) struct PlaneDomainLayout {
    box_xy: [f64; 2],
    mode: PlaneDomainMode,
}

impl PlaneDomainLayout {
    pub(crate) fn new(
        phase: Phase,
        domain_fracs: &[f64],
        box_xy: [f64; 2],
        stripe_axis: usize,
        raft_center: [f64; 2],
        anti_register: bool,
    ) -> Self {
        let active = domain_fracs.iter().filter(|&&f| f > 0.0).count();
        let mode = if active <= 1 || phase == Phase::Random {
            PlaneDomainMode::Random
        } else {
            match phase {
                Phase::Random => PlaneDomainMode::Random,
                Phase::Stripe => {
                    let len = box_xy[stripe_axis];
                    let mut bounds = Vec::with_capacity(domain_fracs.len());
                    let mut cursor = 0.0_f64;
                    for &f in domain_fracs {
                        let next = (cursor + f * len).min(len);
                        bounds.push((cursor, next));
                        cursor = next;
                    }
                    let shift = if anti_register { 0.5 * len } else { 0.0 };
                    PlaneDomainMode::Stripe {
                        axis: stripe_axis,
                        shift,
                        bounds,
                    }
                }
                Phase::Raft => {
                    let f0 = domain_fracs.first().copied().unwrap_or(1.0);
                    let radius = solve_raft_radius(box_xy[0], box_xy[1], f0);
                    let center = if anti_register {
                        [
                            wrap(raft_center[0] + 0.5 * box_xy[0], box_xy[0]),
                            wrap(raft_center[1] + 0.5 * box_xy[1], box_xy[1]),
                        ]
                    } else {
                        raft_center
                    };
                    PlaneDomainMode::Raft { center, radius }
                }
            }
        };
        Self { box_xy, mode }
    }

    pub(crate) fn sample(&self, domain: usize, rng: &mut impl Rng) -> [f64; 2] {
        let [bx, by] = self.box_xy;
        match &self.mode {
            PlaneDomainMode::Random => [rng.random_range(0.0..bx), rng.random_range(0.0..by)],
            PlaneDomainMode::Stripe {
                axis,
                shift,
                bounds,
            } => {
                let len = self.box_xy[*axis];
                let (lo, hi) = bounds.get(domain).copied().unwrap_or((0.0, len));
                let primary = if hi > lo {
                    rng.random_range(lo..hi)
                } else {
                    lo
                };
                let primary = wrap(primary + *shift, len);
                let secondary = rng.random_range(0.0..self.box_xy[1 - *axis]);
                if *axis == 0 {
                    [primary, secondary]
                } else {
                    [secondary, primary]
                }
            }
            PlaneDomainMode::Raft { center, radius } => {
                let half_x = 0.5 * bx;
                let half_y = 0.5 * by;
                let r2 = radius * radius;
                let (dx, dy) = if domain == 0 {
                    if *radius <= half_x.min(half_y) {
                        let r = radius * rng.random::<f64>().sqrt();
                        let theta = rng.random_range(0.0..TAU);
                        let (sin_t, cos_t) = theta.sin_cos();
                        (r * cos_t, r * sin_t)
                    } else {
                        let lx = radius.min(half_x);
                        let ly = radius.min(half_y);
                        let mut picked = (0.0, 0.0);
                        for _ in 0..10_000 {
                            let dx = rng.random_range(-lx..lx);
                            let dy = rng.random_range(-ly..ly);
                            if dx * dx + dy * dy <= r2 {
                                picked = (dx, dy);
                                break;
                            }
                        }
                        picked
                    }
                } else {
                    let mut picked = (half_x, half_y);
                    for _ in 0..100_000 {
                        let dx = rng.random_range(-half_x..half_x);
                        let dy = rng.random_range(-half_y..half_y);
                        if dx * dx + dy * dy > r2 {
                            picked = (dx, dy);
                            break;
                        }
                    }
                    picked
                };
                [wrap(center[0] + dx, bx), wrap(center[1] + dy, by)]
            }
        }
    }
}

/// Pre-resolved lipid template: all force-field lookups done once per species
/// rather than once per placed bead.
struct BeadTemplate<'a> {
    name: &'a str,
    charge: f64,
    mass: f64,
    packing_radius: f64,
    z: f64,
}

/// Place `grid.count` lipids of `composition` in the `box_x` by `box_y` box
/// (uniformly or phase-separated per `options.phase`), append them to
/// `particles` (along with their packing radii and molecule bead ranges), and
/// return the per-component molecule counts.
#[allow(clippy::too_many_arguments)]
fn place_leaflet(
    force_field: &ForceField,
    composition: &Composition,
    box_x: f64,
    box_y: f64,
    grid: LeafletGrid,
    spacing: f64,
    leaflet: Leaflet,
    options: &BuildOptions,
    rng: &mut impl Rng,
    molecule_id: &mut usize,
    atom_id: &mut usize,
    particles: &mut Vec<Particle>,
    radii: &mut Vec<f64>,
    ranges: &mut Vec<std::ops::Range<usize>>,
) -> Result<Vec<usize>, BuildError> {
    let counts = composition.partition(grid.count);
    if grid.count == 0 {
        return Ok(counts);
    }

    let domain_fracs = domain_area_fractions(composition, &counts, force_field);
    let layout = PlaneDomainLayout::new(
        options.phase,
        &domain_fracs,
        [box_x, box_y],
        0,
        [0.5 * box_x, 0.5 * box_y],
        options.anti_register && leaflet == Leaflet::Lower,
    );

    // Resolve each species once: bead properties, packing radii, and chain z offsets.
    let mut templates: Vec<(&str, Vec<BeadTemplate<'_>>, usize, usize)> =
        Vec::with_capacity(counts.len());
    let mut extra_beads = 0usize;
    for (idx, (component, &count)) in composition.components().iter().zip(&counts).enumerate() {
        if count == 0 {
            continue;
        }
        let lipid = force_field
            .lipid(&component.name)
            .ok_or_else(|| BuildError::UnknownLipid(component.name.clone()))?;
        let phosphate = lipid
            .phosphate_index()
            .ok_or_else(|| BuildError::MissingPhosphate(component.name.clone()))?;
        let mut z: f64 = lipid.bonds[..phosphate].iter().map(|b| b.length).sum();
        let last = lipid.beads.len() - 1;
        let mut beads = Vec::with_capacity(lipid.beads.len());
        for (k, bead_name) in lipid.beads.iter().enumerate() {
            let bt = force_field
                .bead_type(bead_name)
                .ok_or_else(|| BuildError::UnknownBead(bead_name.clone()))?;
            beads.push(BeadTemplate {
                name: bt.name.as_str(),
                charge: bt.charge,
                mass: bt.mass,
                packing_radius: POTENTIAL_MINIMUM_FACTOR * bt.interaction.sigma() * 0.5,
                z,
            });
            if k < last {
                z -= lipid.bonds[k].length;
            }
        }
        extra_beads += beads.len() * count;
        templates.push((
            lipid.name.as_str(),
            beads,
            count,
            composition.domain_of(idx),
        ));
    }

    particles.reserve(extra_beads);
    radii.reserve(extra_beads);
    ranges.reserve(grid.count);

    let start = particles.len();
    for (lipid_name, beads, count, domain) in &templates {
        for _ in 0..*count {
            *molecule_id += 1;
            let mol_start = particles.len();
            // Initial lateral placement (uniform or phase-separated); Lloyd's
            // relaxation equalises the density while preserving domain boundaries.
            let [x, y] = layout.sample(*domain, rng);
            for b in beads {
                *atom_id += 1;
                particles.push(Particle {
                    residue: (*lipid_name).to_string(),
                    bead: b.name.to_string(),
                    molecule_id: *molecule_id,
                    atom_id: *atom_id,
                    position: [x, y, b.z],
                    charge: b.charge,
                    mass: b.mass,
                });
                radii.push(b.packing_radius);
            }
            ranges.push(mol_start..particles.len());
        }
    }

    center_leaflet(&mut particles[start..], spacing, leaflet);
    Ok(counts)
}

/// Build a flat bilayer with embedded proteins (or delegate to `Membrane::build`
/// when `proteins` is empty).
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_with_proteins(
    force_field: &ForceField,
    name: &str,
    upper: &Composition,
    lower: &Composition,
    sizing: Sizing,
    options: &BuildOptions,
    proteins: &[crate::embed::ProteinSpec],
    rng: &mut impl Rng,
    frames: Option<&mut Vec<RelaxFrame>>,
) -> Result<crate::embed::EmbeddedSystem, BuildError> {
    if proteins.is_empty() {
        let m = Membrane::build_inner(
            force_field,
            name,
            upper,
            lower,
            sizing,
            options,
            rng,
            frames,
        )?;
        return Ok(crate::embed::EmbeddedSystem::from_membrane(m));
    }

    let upper_spacing = leaflet_spacing(force_field, upper)?;
    let lower_spacing = leaflet_spacing(force_field, lower)?;
    let plan = resolve_geometry(sizing, upper_spacing, lower_spacing)?;

    let mut placed = crate::embed::place_proteins_membrane(proteins, plan.box_x, plan.box_y)?;

    let mut molecule_id = 0;
    let mut atom_id = 0;
    let mut particles = Vec::new();
    let mut radii = Vec::new();
    let mut ranges = Vec::new();

    place_leaflet(
        force_field,
        upper,
        plan.box_x,
        plan.box_y,
        plan.upper,
        upper_spacing,
        Leaflet::Upper,
        options,
        rng,
        &mut molecule_id,
        &mut atom_id,
        &mut particles,
        &mut radii,
        &mut ranges,
    )?;
    let mut upper_lipids = molecule_id;
    place_leaflet(
        force_field,
        lower,
        plan.box_x,
        plan.box_y,
        plan.lower,
        lower_spacing,
        Leaflet::Lower,
        options,
        rng,
        &mut molecule_id,
        &mut atom_id,
        &mut particles,
        &mut radii,
        &mut ranges,
    )?;

    let (bx, by) = (plan.box_x, plan.box_y);
    let carved = crate::embed::carve_lipids(
        &mut particles,
        &mut radii,
        &mut ranges,
        &mut upper_lipids,
        &placed.obstacles,
        0.0,
        |a, b| crate::embed::dist_sq_plane(a, b, bx, by),
    );

    relax_membrane_with_obstacles(
        &mut particles,
        &radii,
        &ranges,
        upper_lipids,
        &placed.obstacles,
        bx,
        by,
        frames,
    );

    let pruned = crate::embed::prune_clashing_lipids(
        &mut particles,
        &mut radii,
        &mut ranges,
        &mut upper_lipids,
        &placed.raw_obstacles,
        |a, b| crate::embed::dist_sq_plane(a, b, bx, by),
    );

    for particle in &mut particles {
        particle.position[0] -= bx / 2.0;
        particle.position[1] -= by / 2.0;
    }
    for prot in &mut placed.instances {
        prot.translate([-bx / 2.0, -by / 2.0, 0.0]);
    }

    let lipid_counts = crate::embed::count_remaining_lipids(&particles, upper, lower);
    let particles = canonicalize(particles, &lipid_counts);

    let mut z_min = f64::INFINITY;
    let mut z_max = f64::NEG_INFINITY;
    for p in &particles {
        z_min = z_min.min(p.position[2]);
        z_max = z_max.max(p.position[2]);
    }
    for prot in &placed.instances {
        for b in &prot.beads {
            z_min = z_min.min(b.position[2]);
            z_max = z_max.max(b.position[2]);
        }
    }
    let z_span = if z_min.is_finite() && z_max.is_finite() {
        z_max - z_min
    } else {
        0.0
    };
    let box_size = [bx, by, z_span + 2.0 * options.padding];

    Ok(crate::embed::assemble_system(
        name,
        placed,
        particles,
        lipid_counts,
        box_size,
        carved,
        pruned,
    ))
}

/// Relax the randomly-placed leaflets to a uniform, clash-free configuration: an
/// independent centroidal-Voronoi (Lloyd) pass per leaflet plus a shared de-clash.
/// `upper_lipids` is the number of lipids in the upper leaflet (placed first), so
/// the two leaflets' molecule ranges can be separated and relaxed at their own
/// densities while the de-clash resolves any overlap between them.
fn relax_membrane(
    particles: &mut [Particle],
    radii: &[f64],
    ranges: &[std::ops::Range<usize>],
    upper_lipids: usize,
    box_x: f64,
    box_y: f64,
    frames: Option<&mut Vec<RelaxFrame>>,
) {
    relax_membrane_with_obstacles(
        particles,
        radii,
        ranges,
        upper_lipids,
        &[],
        box_x,
        box_y,
        frames,
    );
}

#[allow(clippy::too_many_arguments)]
fn relax_membrane_with_obstacles(
    particles: &mut [Particle],
    radii: &[f64],
    ranges: &[std::ops::Range<usize>],
    upper_lipids: usize,
    obstacles: &[pack::Obstacle],
    box_x: f64,
    box_y: f64,
    frames: Option<&mut Vec<RelaxFrame>>,
) {
    if ranges.len() < 2 {
        return;
    }
    let (upper, lower) = ranges.split_at(upper_lipids.min(ranges.len()));
    let leaflets: Vec<&[std::ops::Range<usize>]> = [upper, lower]
        .into_iter()
        .filter(|l| !l.is_empty())
        .collect();
    let mut positions: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();

    // When recording, push one frame per relaxation step. A lipid translates
    // rigidly in x/y, so its first bead tracks its lateral position. The closure
    // owns the per-leaflet ranges and the frame buffer; with no recording the
    // relaxer gets `None` and pays nothing.
    let mut record;
    let observer: Option<pack::StepObserver> = if let Some(frames) = frames {
        let upper = upper.to_vec();
        let lower = lower.to_vec();
        let leaflet_points = |pos: &[[f64; 3]], ranges: &[std::ops::Range<usize>]| {
            ranges
                .iter()
                .map(|r| [pos[r.start][0], pos[r.start][1]])
                .collect()
        };
        record = move |pos: &[[f64; 3]]| {
            frames.push(RelaxFrame {
                upper: leaflet_points(pos, &upper),
                lower: leaflet_points(pos, &lower),
            });
        };
        Some(&mut record)
    } else {
        None
    };

    pack::equalize_plane_with_obstacles(
        &mut positions,
        radii,
        &leaflets,
        obstacles,
        [box_x, box_y],
        DEFAULT_PACK_ITERATIONS,
        DEFAULT_PACK_TOLERANCE,
        observer,
    );
    for (p, pos) in particles.iter_mut().zip(positions) {
        p.position = pos;
    }
}

/// Lift a leaflet so its lowest bead sits half a spacing above the midplane, and
/// mirror it through z when it is the lower leaflet. (Lateral positions are already
/// uniform-random in the box; Lloyd's relaxation arranges the in-plane layout.)
fn center_leaflet(particles: &mut [Particle], spacing: f64, leaflet: Leaflet) {
    let z_min = particles.iter().map(|p| p.position[2]).fold(0.0, f64::min);
    let gap = 0.5 * spacing;
    let z_shift = gap - z_min;
    let factor = match leaflet {
        Leaflet::Upper => 1.0,
        Leaflet::Lower => -1.0,
    };
    for particle in particles {
        particle.position[2] = (particle.position[2] + z_shift) * factor;
    }
}

/// The periodic box (nm): the given `box_x` by `box_y` in plane, and a z extent
/// of the membrane thickness plus `padding` on each side.
fn periodic_box(particles: &[Particle], box_x: f64, box_y: f64, padding: f64) -> [f64; 3] {
    let mut z_min = f64::INFINITY;
    let mut z_max = f64::NEG_INFINITY;
    for particle in particles {
        z_min = z_min.min(particle.position[2]);
        z_max = z_max.max(particle.position[2]);
    }
    [box_x, box_y, (z_max - z_min) + 2.0 * padding]
}

/// Combine per-leaflet counts into total molecule counts per species, ordered by
/// first appearance across the upper then lower compositions.
pub(crate) fn aggregate_counts(
    upper: &Composition,
    upper_counts: &[usize],
    lower: &Composition,
    lower_counts: &[usize],
) -> Vec<(String, usize)> {
    let mut totals: Vec<(String, usize)> = Vec::new();
    for (composition, counts) in [(upper, upper_counts), (lower, lower_counts)] {
        for (component, &count) in composition.components().iter().zip(counts) {
            match totals.iter_mut().find(|(name, _)| *name == component.name) {
                Some(entry) => entry.1 += count,
                None => totals.push((component.name.clone(), count)),
            }
        }
    }
    totals
}

/// Reorder particles so molecules are grouped by species (in `order`) and
/// renumber molecule and atom ids sequentially, moving every particle without
/// cloning. Within a species, molecules keep their placement order, and each
/// molecule's beads stay contiguous and in chain order.
pub(crate) fn canonicalize(
    mut particles: Vec<Particle>,
    order: &[(String, usize)],
) -> Vec<Particle> {
    if order.len() <= 1 {
        let mut molecule_id = 0;
        let mut prev_mol = usize::MAX;
        for (i, p) in particles.iter_mut().enumerate() {
            if p.molecule_id != prev_mol {
                molecule_id += 1;
                prev_mol = p.molecule_id;
            }
            p.molecule_id = molecule_id;
            p.atom_id = i + 1;
        }
        return particles;
    }

    let mut buckets: Vec<Vec<Particle>> = order
        .iter()
        .map(|(_, count)| Vec::with_capacity(count * 6))
        .collect();
    let mut prev_mol = usize::MAX;
    let mut current_bucket: Option<usize> = None;
    for p in particles {
        if p.molecule_id != prev_mol {
            prev_mol = p.molecule_id;
            current_bucket = order.iter().position(|(name, _)| *name == p.residue);
        }
        if let Some(b) = current_bucket {
            buckets[b].push(p);
        }
    }

    let total: usize = buckets.iter().map(Vec::len).sum();
    let mut result = Vec::with_capacity(total);
    let mut molecule_id = 0;
    let mut atom_id = 0;
    for bucket in buckets {
        let mut prev_mol = usize::MAX;
        for mut p in bucket {
            if p.molecule_id != prev_mol {
                molecule_id += 1;
                prev_mol = p.molecule_id;
            }
            atom_id += 1;
            p.molecule_id = molecule_id;
            p.atom_id = atom_id;
            result.push(p);
        }
    }
    result
}

/// Format `particles` into `out`, parallelising across chunks on large systems
/// so coordinate/structure writers saturate all cores without duplicating boilerplate.
pub(crate) fn write_particles(
    out: &mut String,
    particles: &[Particle],
    bytes_per_particle: usize,
    format_one: impl Fn(&mut String, &Particle) + Sync,
) {
    const PAR_CHUNK: usize = 4096;
    if particles.len() >= PAR_CHUNK * 2 {
        let chunks: Vec<String> = particles
            .par_chunks(PAR_CHUNK)
            .map(|chunk| {
                let mut buf = String::with_capacity(chunk.len() * bytes_per_particle);
                for p in chunk {
                    format_one(&mut buf, p);
                }
                buf
            })
            .collect();
        for chunk in chunks {
            out.push_str(&chunk);
        }
    } else {
        for p in particles {
            format_one(out, p);
        }
    }
}

/// Index ranges of `particles` belonging to each molecule. A molecule's beads
/// are stored consecutively, so a change in `molecule_id` starts a new range.
pub(crate) fn molecule_ranges(particles: &[Particle]) -> Vec<std::ops::Range<usize>> {
    let mut ranges = Vec::new();
    let mut start = 0;
    while start < particles.len() {
        let molecule = particles[start].molecule_id;
        let mut end = start + 1;
        while end < particles.len() && particles[end].molecule_id == molecule {
            end += 1;
        }
        ranges.push(start..end);
        start = end;
    }
    ranges
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn force_field() -> ForceField {
        ForceField::isolf()
    }

    fn pure(lipid: &str) -> Composition {
        Composition::from_weights([(lipid.to_string(), 1.0)]).unwrap()
    }

    fn build(upper: &Composition, lower: &Composition, sizing: Sizing, seed: u64) -> Membrane {
        let ff = force_field();
        let mut rng = StdRng::seed_from_u64(seed);
        Membrane::build(
            &ff,
            DEFAULT_NAME,
            upper,
            lower,
            sizing,
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap()
    }

    /// The clash-free minimal-energy leaflet spacing for a pure composition.
    fn spacing(lipid: &str) -> f64 {
        leaflet_spacing(&force_field(), &pure(lipid)).unwrap()
    }

    #[test]
    fn centered_in_box_and_translate_move_the_whole_membrane() {
        let popc = pure("POPC");
        let m = build(
            &popc,
            &popc,
            Sizing::Count {
                upper: 36,
                lower: 36,
            },
            1,
        );
        // `centered_in_box` puts the particle centroid at half the box on every axis.
        let centered = m.centered_in_box();
        let n = centered.particles.len() as f64;
        let centroid = centered.particles.iter().fold([0.0; 3], |mut acc, p| {
            for (a, v) in acc.iter_mut().zip(p.position) {
                *a += v / n;
            }
            acc
        });
        for (c, half) in centroid.into_iter().zip(centered.box_size) {
            assert!((c - half / 2.0).abs() < 1e-9);
        }
        // `translate` is a rigid shift: every particle moves by exactly `shift`,
        // and the box (so the .gro box line) is untouched. The `--center origin`
        // path shifts a centred build by `-box/2` this way.
        let mut shifted = centered.clone();
        let shift = [-centered.box_size[0] / 2.0, 1.5, -2.0];
        shifted.translate(shift);
        for (a, b) in centered.particles.iter().zip(&shifted.particles) {
            for ((&moved, orig), s) in b.position.iter().zip(a.position).zip(shift) {
                assert!((moved - orig - s).abs() < 1e-12);
            }
        }
        assert_eq!(shifted.box_size, centered.box_size);
    }

    #[test]
    fn pure_popc_membrane_has_periodic_box_and_counts() {
        let membrane = build(
            &pure("POPC"),
            &pure("POPC"),
            Sizing::Count {
                upper: 1024,
                lower: 1024,
            },
            1,
        );

        assert_eq!(membrane.total_lipids(), 2048);
        assert_eq!(membrane.lipid_counts, vec![("POPC".to_string(), 2048)]);
        assert_eq!(membrane.particles.len(), 2048 * 6);

        // 1024 is a perfect square ⇒ 32 cells per side of the leaflet spacing.
        assert!((membrane.box_size[0] - 32.0 * spacing("POPC")).abs() < 1e-9);
        assert!((membrane.box_size[1] - 32.0 * spacing("POPC")).abs() < 1e-9);
        assert!(membrane.box_size[2] > 0.0);
    }

    #[test]
    fn spacing_is_the_clash_free_minimal_energy_value() {
        let ff = force_field();
        let popc = pure("POPC");
        let floor = clash_floor(&ff, &popc).unwrap();
        // Just past the bead potential minimum: √1.08 above the clash floor.
        let spacing = leaflet_spacing(&ff, &popc).unwrap();
        assert!((spacing - ENERGY_PACKING_AREA_FACTOR.sqrt() * floor).abs() < 1e-12);
        assert!(spacing > floor);
    }

    #[test]
    fn no_overlaps_under_the_periodic_box() {
        // With the correct periodic box every inter-molecular contact stays near
        // the lipid spacing; the old bounding box collapsed the seam to ~0.
        let membrane = build(
            &pure("POPC"),
            &pure("POPC"),
            Sizing::Count {
                upper: 36,
                lower: 36,
            },
            2,
        );
        let box_size = membrane.box_size;
        let p = &membrane.particles;

        let mut min_r = f64::INFINITY;
        for a in 0..p.len() {
            for b in (a + 1)..p.len() {
                if p[a].molecule_id == p[b].molecule_id {
                    continue;
                }
                let mut d2 = 0.0;
                for ((pa, pb), length) in p[a].position.iter().zip(&p[b].position).zip(&box_size) {
                    let mut d = pa - pb;
                    d -= length * (d / length).round();
                    d2 += d * d;
                }
                min_r = min_r.min(d2.sqrt());
            }
        }
        assert!(
            min_r > 0.5,
            "min image distance {min_r} nm indicates overlap"
        );
    }

    #[test]
    fn asymmetric_leaflet_counts_and_shared_box() {
        // Different counts of the same lipid: the box is sized to the busier
        // (lower) leaflet; the upper sits sparser.
        let membrane = build(
            &pure("POPC"),
            &pure("POPC"),
            Sizing::Count {
                upper: 20,
                lower: 60,
            },
            4,
        );
        assert_eq!(membrane.total_lipids(), 80);
        let (mut upper, mut lower) = (0, 0);
        for range in molecule_ranges(&membrane.particles) {
            if membrane.particles[range.start].position[2] > 0.0 {
                upper += 1;
            } else {
                lower += 1;
            }
        }
        assert_eq!((upper, lower), (20, 60));
        // Box sized to the 60-lipid leaflet: ceil(√60) = 8 cells.
        assert!((membrane.box_size[0] - 8.0 * spacing("POPC")).abs() < 1e-9);
    }

    #[test]
    fn box_fills_with_what_fits() {
        let membrane = build(
            &pure("POPC"),
            &pure("POPC"),
            Sizing::Box { x: 20.0, y: 20.0 },
            11,
        );
        let per_side = (20.0_f64 / spacing("POPC")).floor() as usize;
        assert_eq!(membrane.total_lipids(), 2 * per_side * per_side);
        assert!((membrane.box_size[0] - 20.0).abs() < 1e-9);
    }

    #[test]
    fn rectangular_box_fills_x_and_y_independently() {
        let (x, y) = (24.0, 12.0);
        let membrane = build(&pure("POPC"), &pure("POPC"), Sizing::Box { x, y }, 7);
        let s = spacing("POPC");
        let cells_x = (x / s).floor() as usize;
        let cells_y = (y / s).floor() as usize;
        assert_eq!(membrane.total_lipids(), 2 * cells_x * cells_y);
        assert!((membrane.box_size[0] - x).abs() < 1e-9, "box x");
        assert!((membrane.box_size[1] - y).abs() < 1e-9, "box y");
        // Lipids stay inside the rectangular box.
        for p in &membrane.particles {
            assert!(p.position[0] >= -x / 2.0 - 1.0 && p.position[0] <= x / 2.0 + 1.0);
            assert!(p.position[1] >= -y / 2.0 - 1.0 && p.position[1] <= y / 2.0 + 1.0);
        }
    }

    #[test]
    fn count_in_box_packs_into_the_given_box() {
        let per_side = (10.0_f64 / spacing("POPC")).floor() as usize;
        let n = per_side * per_side / 2;
        let membrane = build(
            &pure("POPC"),
            &pure("POPC"),
            Sizing::CountInBox {
                upper: n,
                lower: n,
                x: 10.0,
                y: 10.0,
            },
            3,
        );
        assert_eq!(membrane.total_lipids(), 2 * n);
        assert!((membrane.box_size[0] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn count_in_box_rejects_more_lipids_than_fit() {
        let ff = force_field();
        let mut rng = StdRng::seed_from_u64(0);
        let leaflet = pure("POPC");
        let error = Membrane::build(
            &ff,
            DEFAULT_NAME,
            &leaflet,
            &leaflet,
            Sizing::CountInBox {
                upper: 100_000,
                lower: 1,
                x: 5.0,
                y: 5.0,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap_err();
        assert!(matches!(error, BuildError::BoxTooSmall { .. }));
    }

    #[test]
    fn unknown_lipid_is_rejected() {
        let ff = force_field();
        let mut rng = StdRng::seed_from_u64(0);
        let leaflet = pure("ZZZZ");
        let error = Membrane::build(
            &ff,
            DEFAULT_NAME,
            &leaflet,
            &leaflet,
            Sizing::Count { upper: 4, lower: 4 },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap_err();
        assert!(matches!(error, BuildError::UnknownLipid(_)));
    }

    #[test]
    fn molecule_and_atom_ids_are_contiguous() {
        let membrane = build(
            &pure("DLPC"),
            &pure("DLPC"),
            Sizing::Count {
                upper: 16,
                lower: 16,
            },
            7,
        );
        for (index, particle) in membrane.particles.iter().enumerate() {
            assert_eq!(particle.atom_id, index + 1);
        }
        assert_eq!(membrane.particles.first().unwrap().molecule_id, 1);
        assert_eq!(
            membrane.particles.last().unwrap().molecule_id,
            membrane.total_lipids()
        );
    }

    #[test]
    fn particles_are_grouped_by_species_in_canonical_order() {
        let upper =
            Composition::from_weights([("DOPC".into(), 1.0), ("DOPS".into(), 1.0)]).unwrap();
        let lower = pure("POPC");
        let membrane = build(
            &upper,
            &lower,
            Sizing::Count {
                upper: 100,
                lower: 100,
            },
            5,
        );

        assert_eq!(
            membrane.lipid_counts,
            vec![
                ("DOPC".to_string(), 50),
                ("DOPS".to_string(), 50),
                ("POPC".to_string(), 100),
            ]
        );

        let residue_blocks: Vec<&str> = membrane
            .particles
            .chunk_by(|a, b| a.molecule_id == b.molecule_id)
            .map(|molecule| molecule[0].residue.as_str())
            .collect();
        let mut expected = Vec::new();
        for (lipid, count) in &membrane.lipid_counts {
            expected.extend(std::iter::repeat_n(lipid.as_str(), *count));
        }
        assert_eq!(residue_blocks, expected);
    }

    #[test]
    fn leaflets_are_mirrored_through_the_midplane() {
        let membrane = build(
            &pure("DPPC"),
            &pure("DPPC"),
            Sizing::Count {
                upper: 36,
                lower: 36,
            },
            3,
        );
        // z has no random component, so the leaflets reach equal and opposite extents.
        let max_z = membrane
            .particles
            .iter()
            .map(|p| p.position[2])
            .fold(f64::MIN, f64::max);
        let min_z = membrane
            .particles
            .iter()
            .map(|p| p.position[2])
            .fold(f64::MAX, f64::min);
        assert!((max_z + min_z).abs() < 1e-9);
        assert!(max_z > 0.0);
    }

    #[test]
    fn z_layout_orders_beads_from_head_to_tail() {
        let membrane = build(
            &pure("DLPC"),
            &pure("DLPC"),
            Sizing::Count { upper: 1, lower: 1 },
            0,
        );
        let upper: Vec<&Particle> = membrane
            .particles
            .iter()
            .filter(|p| p.molecule_id == 1)
            .collect();
        let names: Vec<&str> = upper.iter().map(|p| p.bead.as_str()).collect();
        assert_eq!(names, ["CHO", "PHO", "MID", "DL1", "DL2"]);
        for pair in upper.windows(2) {
            assert!(pair[0].position[2] > pair[1].position[2]);
        }
    }

    #[test]
    fn mixed_composition_partitions_each_leaflet() {
        let upper =
            Composition::from_weights([("DOPC".into(), 1.0), ("DOPS".into(), 1.0)]).unwrap();
        let membrane = build(
            &upper,
            &upper,
            Sizing::Count {
                upper: 100,
                lower: 100,
            },
            5,
        );
        let counts: std::collections::HashMap<_, _> =
            membrane.lipid_counts.iter().cloned().collect();
        assert_eq!(counts["DOPC"], 100);
        assert_eq!(counts["DOPS"], 100);
        assert_eq!(membrane.total_lipids(), 200);
    }

    #[test]
    fn phase_raft_and_stripe_keep_domains_spatially_separated() {
        let ff = force_field();
        let comp = Composition::parse("DPPC=1,DOPC=3").unwrap();

        // Circular raft: DPPC should sit much closer to the centre (0, 0) than DOPC.
        let mut rng = StdRng::seed_from_u64(11);
        let raft = Membrane::build(
            &ff,
            DEFAULT_NAME,
            &comp,
            &comp,
            Sizing::Count {
                upper: 200,
                lower: 200,
            },
            &BuildOptions {
                phase: Phase::Raft,
                ..BuildOptions::default()
            },
            &mut rng,
        )
        .unwrap();
        let mean_r = |res: &str| {
            let (mut s, mut n) = (0.0, 0usize);
            for mol in raft
                .particles
                .chunk_by(|a, b| a.molecule_id == b.molecule_id)
            {
                if mol[0].residue == res {
                    s += mol[0].position[0].hypot(mol[0].position[1]);
                    n += 1;
                }
            }
            s / n as f64
        };
        assert!(
            mean_r("DPPC") < mean_r("DOPC") * 0.65,
            "DPPC raft did not stay central: DPPC r={:.2}, DOPC r={:.2}",
            mean_r("DPPC"),
            mean_r("DOPC")
        );

        // Stripe: contiguous x bands (DPPC in the lower x slab, DOPC in the upper x slab).
        let mut rng = StdRng::seed_from_u64(12);
        let stripe = Membrane::build(
            &ff,
            DEFAULT_NAME,
            &comp,
            &comp,
            Sizing::Count {
                upper: 200,
                lower: 200,
            },
            &BuildOptions {
                phase: Phase::Stripe,
                ..BuildOptions::default()
            },
            &mut rng,
        )
        .unwrap();
        let mean_x = |res: &str| {
            let (mut s, mut n) = (0.0, 0usize);
            for mol in stripe
                .particles
                .chunk_by(|a, b| a.molecule_id == b.molecule_id)
            {
                if mol[0].residue == res {
                    s += mol[0].position[0];
                    n += 1;
                }
            }
            s / n as f64
        };
        assert!(
            mean_x("DOPC") - mean_x("DPPC") > stripe.box_size[0] * 0.35,
            "DPPC/DOPC stripes did not separate: DPPC x={:.2}, DOPC x={:.2}",
            mean_x("DPPC"),
            mean_x("DOPC")
        );

        // Anti-registered raft: upper DPPC near origin, lower DPPC far from origin.
        let mut rng = StdRng::seed_from_u64(13);
        let anti = Membrane::build(
            &ff,
            DEFAULT_NAME,
            &comp,
            &comp,
            Sizing::Count {
                upper: 200,
                lower: 200,
            },
            &BuildOptions {
                phase: Phase::Raft,
                anti_register: true,
                ..BuildOptions::default()
            },
            &mut rng,
        )
        .unwrap();
        let (mut up_r, mut up_n, mut lo_r, mut lo_n) = (0.0, 0usize, 0.0, 0usize);
        for mol in anti
            .particles
            .chunk_by(|a, b| a.molecule_id == b.molecule_id)
        {
            if mol[0].residue != "DPPC" {
                continue;
            }
            let r = mol[0].position[0].hypot(mol[0].position[1]);
            if mol[0].position[2] > 0.0 {
                up_r += r;
                up_n += 1;
            } else {
                lo_r += r;
                lo_n += 1;
            }
        }
        assert!(
            (up_r / up_n as f64) < (lo_r / lo_n as f64) * 0.6,
            "anti-registered raft did not offset leaflets"
        );
    }
}
