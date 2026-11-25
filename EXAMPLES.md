# iSoLF examples

Common builds and the flags behind them. Run `isolf membrane --help`, `isolf vesicle --help`, `isolf tubule --help`, or `isolf protein --help` for every option and default. Each command writes its files into `-o`/`--out`. Add `--top` for the topology, or `--inp` for a full GENESIS setup (which implies `--top`).

- [Lipids and compositions](#lipids-and-compositions)
- [Membranes](#membranes)
- [Vesicles](#vesicles)
- [Tubules](#tubules)
- [Phase-separated domains](#phase-separated-domains)
- [Proteins (standalone and embedded)](#proteins-standalone-and-embedded)
- [Umbrella sampling (amino-acid & peptide/protein PMFs)](#umbrella-sampling-amino-acid--peptideprotein-pmfs)
- [Leaflet balancing](#leaflet-balancing)
- [Centering and padding](#centering-and-padding)
- [GENESIS control files and temperature](#genesis-control-files-and-temperature)
- [Viewing in VMD](#viewing-in-vmd)
- [Animating the relaxation](#animating-the-relaxation)
- [Output files](#output-files)

## Lipids and compositions

A composition lists lipids and how much of each, written `NAME=WEIGHT` and separated by commas. Weights are ratios: plain numbers, fractions, or percentages.

```text
POPC                   only POPC (bare name means weight 1)
POPC=3,POPS=1          3 parts POPC to 1 part POPS  (75% / 25%)
POPC=70%,POPS=30%      the same, as percentages
DPPC=1/2,DOPC=1/2      half and half
```

Repeat a composition flag instead of using commas if you prefer: `-L POPC=3 -L POPS=1`. Use `-L`/`--lipids` when both leaflets (or shells) share the same composition, `--upper`/`--lower` for asymmetric membranes, and `--outer`/`--inner` for asymmetric vesicles.

Run `isolf lipids` at any time to list all 35 supported lipids with their acyl chains, net charges, and bead counts.

## Membranes

A flat bilayer built with `isolf membrane`. Give the lipids and how many per leaflet with `-n`/`--count`:

```bash
isolf membrane --lipids POPC --count 256 -o membrane
```

Use `--upper` and `--lower` for an asymmetric membrane:

```bash
isolf membrane --upper POPC --lower "POPC=3,POPS=1" --count 256 -o asymmetric
```

Set different counts per leaflet with `up=` and `lo=`:

```bash
isolf membrane --upper POPC --lower "POPC=3,POPS=1" --count "up=256,lo=200" -o uneven
```

Or size by in-plane box dimensions in nm with `-b`/`--box` (square or rectangular):

```bash
isolf membrane --lipids POPC --box 20 -o box_20nm
isolf membrane --lipids POPC --box "x=30,y=15" -o box_30x15
```

You can also combine `--count` and `--box` to pack an exact lipid count into a fixed box:

```bash
isolf membrane --lipids POPC --count 200 --box 20 -o sparse_20nm
```

## Vesicles

A closed spherical bilayer built with `isolf vesicle`. Give the composition and `-r`/`--radius` in nm. `ro=` sets the outer radius, `ri=` the inner (lumen) radius, and you can set both. A bare `-r 20` is the same as `-r "ro=20"`.

```bash
isolf vesicle --lipids POPC --radius 20 -o vesicle               # outer radius 20 nm, inner derived
isolf vesicle --lipids POPC --radius "ri=15" -o vesicle_ri       # inner radius 15 nm, outer derived
isolf vesicle --lipids POPC --radius "ri=14,ro=20" -o shell      # both radii fixed
```

Use `--outer` and `--inner` for asymmetric shells:

```bash
isolf vesicle --outer POPC --inner "POPC=3,POPS=1" --radius 20 -o asym_vesicle
```

If you fix both radii and the gap is thinner than the bilayer, `isolf` stops and prints the minimum width it needs.

## Tubules

A cylindrical bilayer periodic along `z` built with `isolf tubule`. Give the composition, `-r`/`--radius` in nm (same `ro=` / `ri=` syntax as `vesicle`), and the periodic axial length `-l`/`--length` in nm:

```bash
isolf tubule --lipids POPC --radius 10 --length 25 -o tubule           # outer radius 10 nm, length 25 nm
isolf tubule --lipids POPC --radius "ri=6" --length 30 -o tubule_ri    # inner lumen radius 6 nm
isolf tubule --outer POPC --inner "POPC=3,POPS=1" --radius "ri=6,ro=10.5" --length 25 -o asym_tubule
```

Like vesicles, tubules generate `nvt.inp` / `pro.inp` under `--inp` so vacuum padding in `x` and `y` is preserved while `z` stays periodic.

## Phase-separated domains

By default lipids are uniformly mixed (`--phase random`). Pass `--phase stripe` or `--phase raft` on `membrane`, `vesicle`, or `tubule` to start with spatially segregated domains at matched local areal density:

- `--phase stripe`: parallel bands along `x` on a flat membrane, latitude bands on a vesicle, or axial rings along `z` on a tubule.
- `--phase raft`: a circular domain (or polar cap on a vesicle) for the first domain surrounded by a mixed matrix of the remaining domains.

With comma-separated species (`"DPPC=1,DOPC=3"`), each species is its own domain. To group multiple lipids into one mixed phase domain, separate domains with `;`:

```bash
# Circular DPPC raft in a DOPC matrix (registered across both leaflets)
isolf membrane -L "DPPC=1,DOPC=3" -n 400 --phase raft --vmd-color-by whole -o raft

# Multi-lipid Lo raft (DPPC+DSPC) in a mixed Ld matrix (POPC+POPS), anti-registered
isolf membrane -L "DPPC=1,DSPC=1;POPC=2,POPS=1" -n 512 --phase raft --anti-register -o antireg_raft

# Latitude-banded vesicle and axially-striped tubule
isolf vesicle -L "DPPC=1,POPC=1" -r 15 --phase stripe --vmd-color-by whole -o banded_vesicle
isolf tubule  -L "DPPC=1,POPC=3" -r 10 -l 25 --phase raft --anti-register --gif -o raft_tubule
```

`--anti-register` offsets the domains on the opposite leaflet/shell (half-period shift for stripes/rings, opposite corner/pole/side for rafts).

## Proteins (standalone and embedded)

### Standalone protein coarse-graining (`isolf protein`)

Convert any all-atom `.pdb` or `.cif`/`.mmcif` protein into a GENESIS RESIDCG AICG2+ (default) or Clementi Gō model. Writes `<name>.gro`, `<name>.itp`, `aicg2p.itp`, `<name>.top`, and `<name>.log` automatically:

```bash
# Default AICG2+ model with GENESIS NVT control files and VMD script
isolf protein 1aki.pdb --inp --vmd -o lysozyme_cg

# Intrinsically disordered N-terminal tail (1-25) + RESPAC charges
isolf protein protein.cif --idr-local "1-25" --idr-nonlocal "1-18" --charges respac.dat -o idr_cg

# Clementi Gō Cα model
isolf protein protein.pdb --model clementi --center origin -o go_cg
```

### Embedding proteins into membranes, vesicles, and tubules (`-p`/`--protein`)

Pass repeatable `-p`/`--protein PATH[,KEY=VAL,...]` on `membrane`, `vesicle`, or `tubule`. By default, input orientation and `z` position are preserved (matching OPM-aligned membrane proteins whose bilayer midplane sits at `z = 0`). Overlapping lipids are carved, remaining lipids relax smoothly around the protein via pinned-site Voronoi CVT and obstacle-aware de-clash, and `--top` / `--inp` emits the combined topology (`aicg2p.itp`, `isolf.itp`, `<prot>.itp`, and `<name>.top`):

```bash
# OPM-oriented channel centered laterally in a 20 x 20 nm bilayer
isolf membrane -L "POPC=3,POPS=1" -b 20 -p channel_opm.pdb --inp --vmd -o channel_mem

# Non-OPM protein centered at the midplane, tilted 15 deg, with 0.2 nm clearance
isolf membrane -L POPC -n 400 -p helix.pdb,z=center,tilt=15,clearance=0.2 --top --vmd -o tilted

# Two different proteins placed side by side
isolf membrane -L POPC -b "x=24,y=14" -p a.pdb,x=-5 -p b.pdb,x=5,rot=90 --top -o dimer_mem

# Multiple copies distributed evenly over a vesicle (Fibonacci sphere) or tubule (helical)
isolf vesicle -L POPC -r 18 -p channel_opm.pdb,copies=6 --inp --vmd -o prot_vesicle
isolf tubule  -L POPC -r 10 -l 30 -p channel_opm.pdb,copies=4 --inp --vmd -o prot_tubule
```

Supported `-p` keys: `z=<nm|center|keep>`, `x=<nm>`, `y=<nm>`, `tilt=<deg>`, `rot=<deg>`, `copies=<N>`, `clearance=<nm>` (default `0.15`), `model=<aicg2+|clementi>`, `scale=<general|average>`, `safe-dihedral=<0|1|3>`, `idr=<ranges>`, `idr-local=<ranges>`, `idr-nonlocal=<ranges>` (separate multiple ranges with `+` or `;`), `charges=<file>`, `name=<molname>`. You can also place a single amino-acid bead directly with `-p aa:TRP,z=1.5[,x=...,y=...,charge=...]`.

## Umbrella sampling (amino-acid & peptide/protein PMFs)

Set up 1D umbrella sampling along the bilayer normal (`z`) for a single amino-acid bead (`TRP`, `TRP,ARG`, `all`) or a multi-residue peptide/protein (`walp23.pdb,tilt=90`) to compute solvation or membrane-insertion PMFs. Each window centers the bilayer at `z = 0` via PBC-safe circular-mean unwrapping, preserves an **identical lipid count across all windows** (zero lipids carved; multi-bead probes use rigid-molecule lateral lipid de-clash), writes a decoupled `lipid_protein.itp` (`[ cg_ISOLF_parameters ]` + Lorentz–Berthelot head `[ cg_WCA_parameters ]` + optional `[ cg_LJ_parameters ]`), GENESIS `DISTMASS` + `direction1 = Z2` control files (`U = k * (dz - z0)^2`), `mbar.inp` for `mbar_analysis`, and `manifest.json`:

```bash
# 41 independent windows from z = 0 to 4.0 nm (step 0.1 nm) for TRP in a 128-lipid POPC bilayer
isolf umbrella setup --probe TRP -L POPC -n 64 --z 0.0..4.0:0.1 -k 0.5 -o us_trp

# Reuse a pre-equilibrated bilayer .gro and run Replica-Exchange Umbrella Sampling (REUS)
isolf umbrella setup --probe TRP --bilayer popc_eq.gro --z 0.0..4.0:0.1 --mode reus \
      --exchange-period 1000 --lipid-protein trial.txt -o reus_trp

# Build campaigns for all 20 canonical amino acids at once (one subfolder per AA)
isolf umbrella setup --probe all --bilayer popc_eq.gro --z 0.0..4.0:0.1 -o campaign_all

# Peptide / transmembrane helix insertion PMF (reusing optimized iSoLFv2 parameters)
isolf umbrella setup --probe walp23.pdb,tilt=90 -L POPC -n 128 \
      --lipid-protein isolfv2_opt.itp --z 0.0..4.5:0.1 --mode reus -o pmf_walp

# Fast Bayesian optimization iteration loop:
# 1) Rewrite only lipid_protein.itp in any umbrella or membrane -p directory (< 1 ms)
isolf params --params trial_02.txt -o us_trp
# 2) Run atdyn across windows (or mpirun -np N atdyn equil_reus.inp / prod_reus.inp)
# 3) Extract RESTR_CVS001 from atdyn logs (and demux REUS .rem files) into cv/{}.dat
isolf umbrella cv --discard 50000 -o us_trp
# 4) Compute PMF with GENESIS mbar_analysis
(cd us_trp && mbar_analysis mbar.inp)
```

Pass `--lipid-protein FILE` on `membrane`, `vesicle`, or `tubule` with `-p` to attach the same `lipid_protein.itp` layer to regular membrane+protein simulations.

Example `--params` / `--lipid-protein` file:

```ini
[cg_ISOLF_parameters]
; name  epsilon(kJ/mol)  sigma(nm)  omega(nm)
TRP     2.200            0.540      1.150

[cg_LJ_parameters]
; optional explicit attraction to lipid head beads (overrides auto WCA for that pair)
TRP PHO 1.800            0.550
```

## Leaflet balancing

In count mode (`--count` without `--box`), if the two leaflets cover different areas, the lighter leaflet is grown to match the denser one so the shared box holds both at matched lateral density. Turn this off with `--no-balance` to keep the exact counts given:

```bash
isolf membrane --upper DOPC --lower DPPC --count "up=256,lo=256" --no-balance -o raw_counts
```

## Centering and padding

By default the finished system sits on the coordinate origin `(0, 0, 0)` (`--center origin`); VMD scripts shift it back into the periodic box for display, so the rendered view is unaffected. Pass `--center box` to centre it in `[0, box]` instead.

`--padding` sets the empty space in nm added on each side (along z for a membrane, along x and y for a tubule, and on every side for a vesicle; 10 nm by default):

```bash
isolf membrane --lipids POPC --count 256 --center box --padding 15 -o padded
```

## GENESIS control files and temperature

`--inp` writes the three-stage GENESIS control files (`min.inp`, `npt.inp` or `nvt.inp`, `pro.inp`) and implies `--top`.

`-t`/`--temperature <K>` sets the thermostat temperature in the control files (323.15 K by default). The starting geometry itself is temperature-independent.

Override per-stage parameters with `--min`, `--equil`, and `--prod` (each accepts a bare step count or comma-separated `KEY=VALUE` pairs, and implies `--inp`):

```bash
isolf membrane --lipids "POPC=4,POPS=1" --count 256 \
      --temperature 310 --equil "steps=50000,dt=0.01" --prod "steps=1000000,crdout=1000" \
      --seed 1 -o body_temp
```

| Stage flag | Keys and defaults |
| --- | --- |
| `--min` | `steps=1000`, `eneout=10`, `crdout=1000`, `rstout=1000`, `nbupdate=10` |
| `--equil` | `steps=20000`, `dt=0.01`, `temp=--temperature`, `eneout=100`, `crdout=100`, `rstout=1000`, `stoptr=10`, `nbupdate=10`, `seed=random` |
| `--prod` | `steps=500000`, `dt=0.02`, `temp=--temperature`, `eneout=500`, `crdout=500`, `rstout=1000`, `stoptr=10`, `nbupdate=10`, `seed=random` |

Placement is random, then relaxed deterministically: pass `--seed` and the same command reproduces byte-identical coordinates. Without `--seed`, `isolf` picks a random one and records it in `<name>.log`.

## Viewing in VMD

`--vmd` writes a `.vmd` script and the `.psf` it loads (passing any `--vmd-*` flag enables `--vmd` automatically). Open it with `vmd -e membrane.vmd`.

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
isolf membrane --lipids "POPC=3,POPS=1" --count 400 \
      --vmd-cutaway --vmd-color-by whole -o viewable
```

With `--inp`, `--vmd` also writes one script per simulation stage (`min.vmd`, `npt.vmd` or `nvt.vmd`, `pro.vmd`) for the trajectories.

## Animating the relaxation

`--gif` writes a `.gif` of the build: Lloyd relaxation spreading the random start, then soft-sphere de-clash settling it (passing any `--gif-*` flag enables `--gif` automatically).

| Option | Default | Effect |
| --- | --- | --- |
| `--gif-mode` | `point` | `point` draws one dot per lipid; `density` draws a local-density heatmap |
| `--gif-scale` | `1` | resolution multiplier (`2` doubles each panel, `0.5` halves it) |
| `--gif-fps` | `14` | playback speed in frames per second |

```bash
isolf membrane --lipids POPC --count 256 --gif -o membrane
isolf vesicle  --lipids POPC --radius 30 --gif-mode density -o vesicle
```

## Output files

Files go in `-o`/`--out`. The base name defaults to the subcommand (`membrane`, `vesicle`, or `tubule`); change it with `--name`.

| File | When | What it is |
| --- | --- | --- |
| `<name>.gro` | always | coordinates in GROMACS format |
| `<name>.log` | always | reproducible run log |
| `<name>.top` | `--top` or `--inp` (always on `protein`) | system topology |
| `isolf.itp` | `--top` or `--inp` (lipid systems) | iSoLF lipid force field |
| `aicg2p.itp`, `<prot>.itp` | `--top`/`--inp` with `-p`, or `protein` | AICG2+ force field and per-protein topology |
| `<name>.pdb` | `--pdb` | coordinates in PDB format |
| `<name>.psf` | `--psf` or `--vmd` | CHARMM structure file |
| `<name>.crd`, `<name>.cif` | `--crd`, `--cif` | CHARMM / mmCIF coordinate formats |
| `min.inp`, `npt.inp`/`nvt.inp`, `pro.inp` | `--inp` | GENESIS control files |
| `<name>.vmd` | `--vmd` | VMD script for the built structure |
| `min.vmd`, … | `--vmd --inp` | VMD script per simulation stage |
| `<name>.gif` | `--gif` | relaxation animation |

A mixed membrane at body temperature, with a fixed seed and all the extra files:

```bash
isolf membrane --lipids "POPC=4,POPS=1" --count 256 \
      --temperature 310 --seed 1 \
      --pdb --inp --vmd -o my_system
```
