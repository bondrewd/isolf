# iSoLF examples

Common builds and the flags behind them. Every option and its default is in `isolf --help`; this page shows how to combine them. Each command writes its files into the folder you give with `--out`. Add `--top` for the topology, or `--inp` for a full GENESIS setup (which implies `--top`).

- [Lipids and compositions](#lipids-and-compositions)
- [Membranes](#membranes)
- [Vesicles](#vesicles)
- [Leaflet balancing](#leaflet-balancing)
- [Centering and padding](#centering-and-padding)
- [Temperature and seed](#temperature-and-seed)
- [Viewing in VMD](#viewing-in-vmd)
- [Animating the relaxation](#animating-the-relaxation)
- [Output files](#output-files)

## Lipids and compositions

A composition lists lipids and how much of each, written `NAME=WEIGHT` and separated by commas. Weights are ratios: plain numbers, fractions, or percentages.

```text
POPC=1                 only POPC
POPC=3,POPS=1          3 parts POPC to 1 part POPS  (75% / 25%)
POPC=70%,POPS=30%      the same, as percentages
DPPC=1/2,DOPC=1/2      half and half
```

A bare name means weight 1, so `--upper POPC` is the same as `--upper POPC=1`. You can also repeat the flag instead of using commas: `--upper POPC=3 --upper POPS=1`.

A lipid name has 4 letters: the first 2 are the tails, the last 2 are the head group.

| Tails | | Head group | |
| --- | --- | --- | --- |
| `DL` | lauroyl (12:0) | `PC` | phosphatidylcholine |
| `DM` | myristoyl (14:0) | `PE` | phosphatidylethanolamine |
| `DP` | palmitoyl (16:0) | `PG` | phosphatidylglycerol |
| `DS` | stearoyl (18:0) | `PS` | phosphatidylserine |
| `DO` | dioleoyl (18:1) | `PA` | phosphatidic acid |
| `PO` | palmitoyl-oleoyl | | |
| `SO` | stearoyl-oleoyl | | |

So `POPC`, `DOPS`, `DPPE`, and the rest: 35 lipids in all. Common choices are POPC, DOPC, DPPC, POPE, and POPS.

## Membranes

A flat bilayer. Give the lipids and how many per leaflet:

```bash
isolf --upper POPC=1 --lipids-per-leaflet 256 --out membrane
```

Add `--lower` for a different bottom leaflet. Cell membranes are usually asymmetric:

```bash
isolf --upper POPC=1 --lower "POPC=3,POPS=1" --lipids-per-leaflet 256 --out asymmetric
```

Set the two counts apart with `up=` and `lo=`:

```bash
isolf --upper POPC=1 --lower "POPC=3,POPS=1" --lipids-per-leaflet "up=256,lo=200" --out uneven
```

Or give the box size in nm and let `isolf` fill it, as a square or a rectangle:

```bash
isolf --upper POPC=1 --membrane 20 --out box_20nm
isolf --upper POPC=1 --membrane "x=30,y=15" --out box_30x15
```

## Vesicles

A closed sphere. Give the composition and a radius in nm. `ro=` sets the outer radius, `ri=` the inner (lumen) radius, and you can give both. A bare `--vesicle 20` is the same as `--vesicle "ro=20"`.

```bash
isolf --upper POPC=1 --vesicle 20 --out vesicle           # outer radius 20 nm, inner found for you
isolf --upper POPC=1 --vesicle "ri=15" --out vesicle_ri   # inner radius 15 nm, outer found for you
isolf --upper POPC=1 --vesicle "ri=14,ro=20" --out shell  # both radii fixed
```

`--upper` is the outer shell and `--lower` the inner (lumen-facing) shell. If you fix both radii and the gap is thinner than the bilayer, `isolf` stops and prints the minimum width it needs.

## Leaflet balancing

When you set both leaflet counts and the two leaflets cover different areas, the lighter leaflet is grown to match the denser one, so the shared box holds both. This is on by default. Turn it off with `--no-balance` to keep the counts exactly as given.

```bash
isolf --upper DOPC=1 --lower DPPC=1 --lipids-per-leaflet "up=256,lo=256" --no-balance --out raw_counts
```

## Centering and padding

By default the finished system sits on the coordinate origin `(0, 0, 0)` (`--center origin`); the VMD scripts shift it back into the periodic box for display, so the rendered view is unaffected. Pass `--center box` to center it in the box `[0, box]` instead. Either way only the coordinates move, the box is unchanged, and the two modes differ by exactly half the box.

`--padding` sets the empty space in nm added along z for a membrane, or on every side for a vesicle (10 nm by default).

```bash
isolf --upper POPC=1 --lipids-per-leaflet 256 --center box --padding 15 --out padded
```

## Temperature and seed

`--temperature <K>` sets the thermostat temperature written into the GENESIS control files (323.15 K by default). The geometry is independent of it: lipids start at a fixed clash-free spacing that an NPT run compresses to the real density.

Placement is random, then relaxed to an even, clash-free layout. The build is reproducible: pass `--seed` and the same command gives byte-identical coordinates. Without `--seed`, `isolf` picks a random one and records it in the log, so you can reproduce any run from there.

```bash
isolf --upper POPC=1 --lipids-per-leaflet 256 --temperature 310 --seed 1 --out body_temp
```

A per-stage flag overrides `--temperature` for one control file: `--nvt-temperature`, `--npt-temperature`, `--pro-temperature`. The step counts and output periods for each stage are in `isolf --help`.

## Viewing in VMD

`--vmd` writes a `.vmd` script and the `.psf` it loads. Open it with `vmd -e membrane.vmd`. A few flags change the look:

| Option | Default | Effect |
| --- | --- | --- |
| `--vmd-color-by` | `head` | `head` tints the head bead by species, `bead` colours every bead by role, `whole` colours the whole lipid by species |
| `--vmd-surface` | off | a smooth QuickSurf surface instead of VDW spheres |
| `--vmd-cutaway` | off | a clipping plane that reveals the cross-section |
| `--vmd-render` | off | append a Tachyon ray-tracing command to each script |
| `--vmd-background` | `white` | a VMD colour name |
| `--vmd-material` | `AOChalky` | a VMD material name, e.g. `Goodsell` |
| `--vmd-smooth` | `5` | trajectory smoothing window, in frames |

```bash
isolf --upper "POPC=3,POPS=1" --lipids-per-leaflet 400 \
      --vmd --vmd-cutaway --vmd-color-by whole --out viewable
```

With `--inp`, `--vmd` also writes one script per simulation step (`min.vmd`, `npt.vmd` or `nvt.vmd`, `pro.vmd`) for the trajectories.

## Animating the relaxation

`--gif` writes a `.gif` of the build: the Lloyd relaxation spreading the random start to an even layout, then the soft-sphere de-clash settling it. The upper leaflet (the outer shell of a vesicle) is the left panel, the lower leaflet (the inner shell) the right.

```bash
isolf --upper POPC=1 --lipids-per-leaflet 256 --gif --out membrane
```

Three flags change it:

| Option | Default | Effect |
| --- | --- | --- |
| `--gif-mode` | `point` | `point` draws one dot per lipid; `density` draws a heatmap of the local density |
| `--gif-scale` | `1` | resolution multiplier, so `2` doubles each panel and `0.5` halves it |
| `--gif-fps` | `14` | playback speed in frames per second |

`point` reads clearly up to a few thousand lipids per leaflet. Past that the dots overlap into a filled disc, so use `density`: it colours each pixel by how many lipids sit near it and stays legible at any count. The heatmap fills every pixel, so its file is several times larger than the point scatter, and `--gif-scale 0.5` brings it back down.

```bash
isolf --upper POPC=1 --vesicle 30 --gif --gif-mode density --out vesicle
```

A vesicle is shown as a disc: each shell is flattened with an equal-area projection, so an even spread of lipids on the sphere stays even on the disc.

## Output files

Files go in the `--out` folder. The base name is the build mode (`membrane` or `vesicle`); change it with `--name`.

| File | When | What it is |
| --- | --- | --- |
| `<name>.gro` | always | the structure (coordinates) |
| `<name>.log` | always | a run log for reproducibility (see below) |
| `<name>.top` | with `--top` or `--inp` | the topology |
| `isolf.itp` | with `--top` or `--inp` | the iSoLF lipid force field |
| `<name>.pdb` | with `--pdb` | coordinates in PDB format |
| `<name>.psf` | with `--psf` (or `--vmd`) | a CHARMM structure file |
| `<name>.crd`, `<name>.cif` | with `--crd`, `--cif` | other coordinate formats |
| `min.inp`, `npt.inp`/`nvt.inp`, `pro.inp` | with `--inp` | GENESIS control files |
| `<name>.vmd` | with `--vmd` | a VMD script for the built structure |
| `<name>.gif` | with `--gif` | an animation of the leaflet relaxation, Lloyd then de-clash |
| `min.vmd`, … | with `--vmd --inp` | a VMD script per simulation step |

`--inp` writes the control files for a 3-step run (minimization, equilibration, production) and implies `--top`, since a run reads the topology. A membrane equilibrates under NPT (`npt.inp`); a vesicle under NVT (`nvt.inp`).

`<name>.log` is written on every run. It records the exact command with the seed pinned in, plus the time, host, user, working directory, and what was built. To reproduce a build, copy the `command` line from the log and run it again.

A mixed membrane at body temperature, with a fixed seed and all the extra files:

```bash
isolf --upper "POPC=4,POPS=1" --lipids-per-leaflet 256 \
      --temperature 310 --seed 1 \
      --pdb --inp --vmd --out my_system
```
