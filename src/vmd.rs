//! VMD (1.9.x) visualization scripts for a built system.
//!
//! [`vmd_script`] renders a script from a [`VmdSource`]: the initial built
//! structure (`<base>.gro`) or a simulation stage's `.dcd` trajectory. Either way
//! it loads the system's `.psf` and draws a representation per bead role *present
//! in the system* — at physical bead radii taken from the iSoLF σ — with ambient-
//! occlusion lighting. The optional touches (a smooth QuickSurf surface, a
//! cut-away clipping plane, a Tachyon render) and the cosmetic knobs (background,
//! material, smoothing window) are driven by [`VmdOptions`], i.e. the `--vmd-*`
//! command-line flags.

use std::fmt::Write;

use crate::force_field::{ForceField, Interaction};
use crate::membrane::Particle;

/// nm → ångström, with the factor of a half that turns a σ into a bead radius.
const SIGMA_TO_RADIUS: f64 = 5.0;

/// The system shape, which selects the camera and PBC handling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    /// Flat bilayer: a side view of the leaflets, with the periodic box drawn.
    Membrane,
    /// Vesicle: a perspective view centred on the sphere.
    Vesicle,
}

/// How lipid beads are coloured in the VMD script (the `--vmd-color-by` scope).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LipidColoring {
    /// Colour each bead by its role (head/phosphate/glycerol/tail).
    Role,
    /// Colour the head bead by lipid species (`ResName`); the rest by role.
    Head,
    /// Colour every lipid bead by lipid species (`ResName`).
    Whole,
}

/// Visualization choices shared by every stage's script (the `--vmd-*` flags).
#[derive(Debug, Clone, Copy)]
pub struct VmdOptions<'a> {
    /// Use a smooth QuickSurf surface instead of VDW spheres.
    pub surface: bool,
    /// Slice every representation with a clipping plane (reveals the cross-section).
    pub cutaway: bool,
    /// Append a Tachyon ray-tracing render command.
    pub render: bool,
    /// Display background colour (a VMD colour name).
    pub background: &'a str,
    /// Representation material (a VMD material name).
    pub material: &'a str,
    /// Trajectory smoothing window, in frames.
    pub smooth: usize,
    /// How lipids are coloured: by bead role (the default), or by lipid species
    /// (`ResName`) for just the head or the whole lipid, so a mixture of lipids is
    /// distinguishable.
    pub lipid_coloring: LipidColoring,
    /// Translate the loaded coordinates into the periodic box before drawing it
    /// (for `--center origin`, whose files are centred on the origin and would
    /// otherwise draw in a corner). The shift comes from VMD's own loaded cell, so
    /// it is unit-correct whatever scaling the coordinate reader applied.
    pub recenter_to_box: bool,
}

/// What a script visualizes.
#[derive(Debug, Clone, Copy)]
pub enum VmdSource<'a> {
    /// The initial built coordinates from `<base>.gro` (a single frame).
    Structure,
    /// Simulation stage `name`'s `<name>.dcd` trajectory, with `frames` frames
    /// (its `nsteps / crdout_period`) when known.
    Trajectory {
        name: &'a str,
        frames: Option<usize>,
    },
}

/// The coarse role a bead plays, which sets its representation colour.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Role {
    Head,
    Phosphate,
    Glycerol,
    Tail,
}

impl Role {
    /// Roles in display order.
    const ORDER: [Role; 4] = [Role::Head, Role::Phosphate, Role::Glycerol, Role::Tail];

    fn label(self) -> &'static str {
        match self {
            Role::Head => "head",
            Role::Phosphate => "phosphate",
            Role::Glycerol => "glycerol",
            Role::Tail => "tail",
        }
    }

    /// VMD `ColorID` for this role.
    fn color_id(self) -> u8 {
        match self {
            Role::Head => 0,      // blue
            Role::Phosphate => 1, // red
            Role::Glycerol => 3,  // orange
            Role::Tail => 6,      // silver
        }
    }
}

/// Classify a bead: phosphate and glycerol are named explicitly; acyl-tail beads
/// carry the iSoLF potential; the rest are head beads.
fn role_of(name: &str, force_field: &ForceField) -> Role {
    match name {
        "PHO" | "PHA" => Role::Phosphate,
        "MID" => Role::Glycerol,
        other => match force_field.bead_type(other) {
            Some(bead) if matches!(bead.interaction, Interaction::Isolf { .. }) => Role::Tail,
            _ => Role::Head,
        },
    }
}

/// Render a VMD script for `source`: it loads `<base>.psf`, appends the initial
/// `<base>.gro` or the stage's `<name>.dcd`, and (for a trajectory) wraps its
/// frames.
pub fn vmd_script(
    force_field: &ForceField,
    particles: &[Particle],
    base: &str,
    layout: Layout,
    source: VmdSource,
    options: &VmdOptions,
) -> String {
    // What to load, what to title it, and the render filename, by source.
    let (label, coord_file, render_name, frames) = match source {
        VmdSource::Structure => (
            "structure".to_string(),
            format!("mol addfile ./{base}.gro type gro waitfor all"),
            base.to_string(),
            None,
        ),
        VmdSource::Trajectory { name, frames } => (
            format!("{name} stage"),
            format!("mol addfile ./{name}.dcd type dcd waitfor all"),
            format!("{base}_{name}"),
            frames,
        ),
    };
    let is_trajectory = matches!(source, VmdSource::Trajectory { .. });

    // Bead names actually present in the system, in a stable sorted order.
    let mut present: Vec<&str> = particles.iter().map(|p| p.bead.as_str()).collect();
    present.sort_unstable();
    present.dedup();

    // Present beads grouped into the roles that have at least one bead.
    let groups: Vec<(Role, Vec<&str>)> = Role::ORDER
        .into_iter()
        .filter_map(|role| {
            let beads: Vec<&str> = present
                .iter()
                .copied()
                .filter(|&name| role_of(name, force_field) == role)
                .collect();
            (!beads.is_empty()).then_some((role, beads))
        })
        .collect();

    let mut s = String::new();
    let shape = match layout {
        Layout::Membrane => "membrane",
        Layout::Vesicle => "vesicle",
    };
    let _ = writeln!(s, "# isolf visualization — {label} ({shape})");
    let _ = writeln!(s, "mol new ./{base}.psf type psf waitfor all");
    let _ = writeln!(s, "{coord_file}");
    let _ = writeln!(s, "mol delrep 0 top\n");

    // `--center origin`: the coordinates are centred on the origin, so shift them
    // into the drawn periodic box (every frame). The half-box comes from VMD's
    // loaded cell, matching the coordinates' units exactly; this makes the view
    // identical to a box-centred build, after which `pbc box`/`pbc wrap` apply as
    // usual.
    if options.recenter_to_box {
        let _ = writeln!(s, "# --center origin: shift the system into the box.");
        let _ = writeln!(s, "set isolf_cell [molinfo top get {{a b c}}]");
        let _ = writeln!(
            s,
            "set isolf_half [vecscale 0.5 $isolf_cell]\nset isolf_sel [atomselect top all]"
        );
        let _ = writeln!(
            s,
            "for {{set isolf_f 0}} {{$isolf_f < [molinfo top get numframes]}} {{incr isolf_f}} {{"
        );
        let _ = writeln!(s, "  $isolf_sel frame $isolf_f");
        let _ = writeln!(s, "  $isolf_sel moveby $isolf_half");
        let _ = writeln!(s, "}}\n$isolf_sel delete\n");
    }

    // Physical bead radii from σ (the .psf carries no radii for CG bead names).
    let _ = writeln!(s, "# Physical bead radii from the iSoLF sigma (angstrom).");
    for &name in &present {
        let radius = if let Some(bead) = force_field.bead_type(name) {
            bead.interaction.sigma() * SIGMA_TO_RADIUS
        } else {
            continue;
        };
        let _ = writeln!(
            s,
            "set sel [atomselect top \"name {name}\"]; $sel set radius {radius:.2}; $sel delete"
        );
    }

    // Representation style: spheres by default, a smooth surface with --vmd-surface.
    let repstyle = if options.surface {
        "QuickSurf 1.0 0.5 1.0 1.0"
    } else {
        "VDW 1.2 10.0"
    };
    let _ = writeln!(s, "\nset repstyle {{{repstyle}}}");

    // One representation per present bead role.
    let _ = writeln!(
        s,
        "\n# One representation per bead role present in the system."
    );
    for (role, beads) in &groups {
        let selection = beads
            .iter()
            .map(|b| format!("name {b}"))
            .collect::<Vec<_>>()
            .join(" or ");
        let _ = writeln!(s, "# {}", role.label());
        let _ = writeln!(s, "mol representation {{*}}$repstyle");
        // Lipid beads can be coloured per lipid type so a mixture is
        // distinguishable: the whole lipid (`Whole`) or just the head (`Head`).
        let by_lipid_type = match options.lipid_coloring {
            LipidColoring::Role => false,
            LipidColoring::Head => *role == Role::Head,
            LipidColoring::Whole => true,
        };
        if by_lipid_type {
            let _ = writeln!(s, "mol color ResName");
        } else {
            let _ = writeln!(s, "mol color ColorID {}", role.color_id());
        }
        let _ = writeln!(s, "mol selection {{{selection}}}");
        let _ = writeln!(s, "mol material {}", options.material);
        let _ = writeln!(s, "mol addrep top");
    }

    // Periodic box (drawn for every source). A trajectory is also wrapped over
    // its frame range; the initial structure is a single, already-centred frame.
    let _ = writeln!(
        s,
        "\n# Periodic box{}.",
        if is_trajectory { " and wrapping" } else { "" }
    );
    let _ = writeln!(s, "pbc box");
    if is_trajectory {
        match frames {
            Some(frames) if frames >= 1 => {
                let _ = writeln!(s, "# {frames} frames = nsteps / crdout_period");
                let _ = writeln!(s, "pbc wrap -center com -first 0 -last {}", frames - 1);
            }
            _ => {
                let _ = writeln!(s, "pbc wrap -center com -all");
            }
        }
    }

    // Optional cut-away clipping plane through the middle of every representation.
    if options.cutaway {
        let _ = writeln!(s, "\n# Cut-away clipping plane (--vmd-cutaway).");
        let _ = writeln!(
            s,
            "for {{set r 0}} {{$r < [molinfo top get numreps]}} {{incr r}} {{"
        );
        let _ = writeln!(s, "  mol clipplane center 0 $r top {{0 0 0}}");
        let _ = writeln!(s, "  mol clipplane normal 0 $r top {{0 1 0}}");
        let _ = writeln!(s, "  mol clipplane status 0 $r top 2");
        let _ = writeln!(s, "}}");
    }

    // Camera.
    let projection = match layout {
        Layout::Membrane => "orthographic",
        Layout::Vesicle => "perspective",
    };
    let _ = writeln!(s, "\n# Camera.");
    if layout == Layout::Membrane {
        let _ = writeln!(s, "rotate x by -90");
    }
    let _ = writeln!(s, "scale by 1.2");

    // Lighting and background.
    let _ = writeln!(s, "\n# Lighting and background.");
    let _ = writeln!(s, "display projection {projection}");
    let _ = writeln!(s, "display depthcue off");
    let _ = writeln!(s, "display shadows on");
    let _ = writeln!(s, "display ambientocclusion on");
    let _ = writeln!(s, "display dof on");
    let _ = writeln!(s, "color Display Background {}", options.background);
    let _ = writeln!(s, "axes location off");

    // Smooth every representation's trajectory (only meaningful with frames).
    if is_trajectory {
        let _ = writeln!(s, "\n# Smooth the trajectory.");
        for rep in 0..groups.len() {
            let _ = writeln!(s, "mol smoothrep 0 {rep} {}", options.smooth);
        }
    }

    let _ = writeln!(s, "\ndisplay resetview");

    // Optional ray-traced still.
    if options.render {
        let _ = writeln!(s, "\n# Ray-traced still (--vmd-render).");
        let _ = writeln!(s, "render TachyonInternal {render_name}.tga");
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composition::Composition;
    use crate::membrane::{BuildOptions, DEFAULT_NAME, Membrane, Sizing};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn popc_membrane() -> Membrane {
        let ff = ForceField::isolf();
        let leaflet = Composition::parse("POPC=1").unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        Membrane::build(
            &ff,
            DEFAULT_NAME,
            &leaflet,
            &leaflet,
            Sizing::Count {
                upper: 16,
                lower: 16,
            },
            &BuildOptions::default(),
            &mut rng,
        )
        .unwrap()
    }

    fn options() -> VmdOptions<'static> {
        VmdOptions {
            surface: false,
            cutaway: false,
            render: false,
            background: "white",
            material: "AOChalky",
            smooth: 5,
            lipid_coloring: LipidColoring::Role,
            recenter_to_box: false,
        }
    }

    fn script(frames: Option<usize>, options: &VmdOptions) -> String {
        vmd_script(
            &ForceField::isolf(),
            &popc_membrane().particles,
            "membrane",
            Layout::Membrane,
            VmdSource::Trajectory {
                name: "npt",
                frames,
            },
            options,
        )
    }

    #[test]
    fn structure_source_loads_the_gro_without_a_trajectory() {
        let script = vmd_script(
            &ForceField::isolf(),
            &popc_membrane().particles,
            "membrane",
            Layout::Membrane,
            VmdSource::Structure,
            &options(),
        );
        assert!(script.contains("mol new ./membrane.psf type psf"));
        assert!(script.contains("mol addfile ./membrane.gro type gro"));
        assert!(!script.contains(".dcd"));
        assert!(!script.contains("pbc wrap")); // a single, already-centred frame
        assert!(!script.contains("mol smoothrep")); // no trajectory to smooth
        assert!(script.contains("pbc box")); // the box is still drawn
    }

    #[test]
    fn loads_psf_dcd_and_wraps_the_computed_frame_range() {
        let script = script(Some(100), &options());
        assert!(script.contains("mol new ./membrane.psf type psf"));
        assert!(script.contains("mol addfile ./npt.dcd type dcd"));
        assert!(script.contains("color Display Background white"));
        assert!(script.contains("mol material AOChalky"));
        // 100 frames ⇒ indices 0..99.
        assert!(script.contains("# 100 frames = nsteps / crdout_period"));
        assert!(script.contains("pbc wrap -center com -first 0 -last 99"));
    }

    #[test]
    fn defaults_omit_the_optional_visualizations() {
        let script = script(Some(100), &options());
        assert!(script.contains("set repstyle {VDW 1.2 10.0}"));
        assert!(!script.contains("QuickSurf"));
        assert!(!script.contains("clipplane"));
        assert!(!script.contains("render TachyonInternal"));
    }

    #[test]
    fn flags_enable_surface_cutaway_render_and_set_cosmetics() {
        let opts = VmdOptions {
            surface: true,
            cutaway: true,
            render: true,
            background: "black",
            material: "Goodsell",
            smooth: 9,
            lipid_coloring: LipidColoring::Role,
            recenter_to_box: false,
        };
        let script = script(Some(40), &opts);
        assert!(script.contains("set repstyle {QuickSurf 1.0 0.5 1.0 1.0}"));
        assert!(script.contains("mol clipplane status 0 $r top 2"));
        assert!(script.contains("render TachyonInternal membrane_npt.tga"));
        assert!(script.contains("color Display Background black"));
        assert!(script.contains("mol material Goodsell"));
        assert!(script.contains("mol smoothrep 0 0 9"));
    }

    #[test]
    fn lipid_coloring_head_colours_only_heads_by_resname() {
        let plain = script(Some(50), &options());
        assert!(plain.contains("mol color ColorID 0")); // head role is colour-id by default
        let by_lipid = script(
            Some(50),
            &VmdOptions {
                lipid_coloring: LipidColoring::Head,
                ..options()
            },
        );
        assert!(by_lipid.contains("mol color ResName"));
        // Tails keep their role colour even when heads are per-lipid.
        assert!(by_lipid.contains("mol color ColorID 6"));
    }

    #[test]
    fn lipid_coloring_whole_colours_every_lipid_bead_by_resname() {
        // The whole lipid by type: even the tail (ColorID 6 by default) becomes
        // ResName, so a membrane-only script keeps no role colour at all.
        let plain = script(Some(50), &options());
        assert!(plain.contains("mol color ColorID 6")); // tail is colour-id by default
        let by_lipid = script(
            Some(50),
            &VmdOptions {
                lipid_coloring: LipidColoring::Whole,
                ..options()
            },
        );
        assert!(by_lipid.contains("mol color ResName"));
        assert!(!by_lipid.contains("mol color ColorID")); // no lipid role keeps a role colour
    }

    #[test]
    fn recenter_to_box_shifts_the_coordinates_into_the_box() {
        // --center origin: the structure script (which never wraps) must translate
        // the origin-centred coordinates into the drawn box, sized from VMD's own
        // loaded cell, and do so before `pbc box` frames the system.
        let plain = vmd_script(
            &ForceField::isolf(),
            &popc_membrane().particles,
            "membrane",
            Layout::Membrane,
            VmdSource::Structure,
            &options(),
        );
        assert!(!plain.contains("moveby")); // box mode shifts nothing
        let shifted = vmd_script(
            &ForceField::isolf(),
            &popc_membrane().particles,
            "membrane",
            Layout::Membrane,
            VmdSource::Structure,
            &VmdOptions {
                recenter_to_box: true,
                ..options()
            },
        );
        assert!(shifted.contains("molinfo top get {a b c}"));
        assert!(shifted.contains("$isolf_sel moveby $isolf_half"));
        assert!(shifted.find("moveby").unwrap() < shifted.find("pbc box").unwrap());
    }

    #[test]
    fn only_present_beads_get_representations() {
        let script = script(Some(50), &options());
        for present in ["name CHO", "name PHO", "name MID", "name PO1"] {
            assert!(script.contains(present), "missing {present}");
        }
        for absent in ["name GLC", "name SRI", "name ETH", "name DL1", "name DO1"] {
            assert!(!script.contains(absent), "unexpected {absent}");
        }
    }

    #[test]
    fn no_coordinate_output_falls_back_to_wrapping_all_frames() {
        let script = script(None, &options());
        assert!(script.contains("pbc wrap -center com -all"));
    }

    #[test]
    fn vesicle_uses_a_perspective_view() {
        let script = vmd_script(
            &ForceField::isolf(),
            &popc_membrane().particles,
            "vesicle",
            Layout::Vesicle,
            VmdSource::Trajectory {
                name: "pro",
                frames: Some(500),
            },
            &options(),
        );
        assert!(script.contains("display projection perspective"));
        assert!(script.contains("pbc box")); // the periodic box is drawn for every geometry
        assert!(!script.contains("rotate x by -90"));
    }
}
