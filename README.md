<div align="center">

<img src="assets/banner.svg" alt="iSoLF" width="660">

**Coarse-grained starting structures for [GENESIS](https://mdgenesis.org/) molecular dynamics: lipid membranes, vesicles, and tubules.**

[![release](https://img.shields.io/github/v/release/bondrewd/isolf?sort=semver&color=2dd4bf&label=release)](https://github.com/bondrewd/isolf/releases)
[![ci](https://img.shields.io/github/actions/workflow/status/bondrewd/isolf/ci.yml?branch=main&color=3b82f6&label=ci)](https://github.com/bondrewd/isolf/actions)
[![downloads](https://img.shields.io/github/downloads/bondrewd/isolf/total?color=8b5cf6&label=downloads)](https://github.com/bondrewd/isolf/releases)
[![license](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-64748b)](#-license)

</div>

You give `isolf` a lipid recipe, and it writes the coordinates. With `--top` or `--inp` it also writes the topology and force-field files GENESIS reads. Lipids use the iSoLF implicit-solvent model, so there's no water to add. It's a single binary, with nothing else to install.

<p align="center"><img src="assets/demo.gif" alt="iSoLF building a membrane" width="760"></p>

## 🗺️ Map

[Install](#-install) · [Quick start](#-quick-start) · [Commands](#-commands) · [Examples](EXAMPLES.md) · [Output](#-output) · [Run in GENESIS](#-run-in-genesis) · [License](#-license)

## 📦 Install

**macOS / Linux**

```bash
curl -LsSf https://raw.githubusercontent.com/bondrewd/isolf/main/install.sh | sh
```

**Windows (PowerShell)**

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://raw.githubusercontent.com/bondrewd/isolf/main/install.ps1 | iex"
```

The installer picks the right binary for your machine, checks its sha256, drops it in `~/.local/bin`, and puts that on your PATH. Update in place with `isolf update`. Remove it with `isolf uninstall`.

<details>
<summary>Other ways to install</summary>

```bash
# with wget instead of curl
wget -qO- https://raw.githubusercontent.com/bondrewd/isolf/main/install.sh | sh

# pin a version, or change the install directory
curl -LsSf https://raw.githubusercontent.com/bondrewd/isolf/main/install.sh | sh -s -- 0.1.0
curl -LsSf https://raw.githubusercontent.com/bondrewd/isolf/main/install.sh | ISOLF_INSTALL_DIR=/usr/local/bin sh

# from source, with a Rust toolchain (https://rustup.rs)
cargo install --path .
```

You can also download a prebuilt binary from the [Releases page](https://github.com/bondrewd/isolf/releases). Pick `windows-msvc.zip`, `apple-darwin.tar.gz`, or one of the `*-linux-musl.tar.gz` files, unzip it, and run `./isolf --version`. On macOS, clear the quarantine flag once with `xattr -dr com.apple.quarantine ./isolf`.

</details>

## 🚀 Quick start

```bash
isolf membrane --lipids POPC --count 256 --top -o my_first_membrane
```

This writes `my_first_membrane/` with the structure (`membrane.gro`), the topology (`membrane.top`), and the lipid force field (`isolf.itp`). Add `--inp` for a full GENESIS setup, or `--vmd` for a viewer script.

## 🧱 Commands

| Command | Example |
| --- | --- |
| **`membrane`** | `isolf membrane -L "POPC=3,POPS=1" -n 256 -o m` |
| **`vesicle`** | `isolf vesicle -L POPC -r 20 -o v` |
| **`tubule`** | `isolf tubule -L POPC -r 10 -l 25 -o t` |
| **`lipids`** | `isolf lipids` (lists all 35 supported lipids) |

Lipids are `NAME=WEIGHT` recipes like `"POPC=3,POPS=1"` (repeat or comma-separate to combine). Use `-L`/`--lipids` for symmetric systems, `--upper`/`--lower` for asymmetric membranes, and `--outer`/`--inner` for asymmetric vesicles and tubules. Add `--phase stripe` or `--phase raft` (with `;`-separated multi-lipid domains and optional `--anti-register`) to build phase-separated domains on any geometry. Size a membrane by lipid count (`-n`/`--count`), box size (`-b`/`--box`, e.g. `x=10,y=20`), or both; size a vesicle by `-r`/`--radius` (`20`, `ro=20`, `ri=15`, or `ri=14,ro=20`); size a tubule by `-r`/`--radius` and periodic `-l`/`--length`.

While it runs, `isolf` prints each phase and a closing summary: the composition, the box, the particle count, the time, and the command to run next. Use `-q` for the summary alone, `-v` for full file paths, `--ascii` for plain glyphs.

> 🧪 More recipes are in **[EXAMPLES.md](EXAMPLES.md)**: asymmetric leaflets, rectangular membranes, vesicles, tubules, phase-separated domains, MD stage overrides, and output formats.

## 📦 Output

Everything lands in the `-o`/`--out` folder. The `.gro` structure and a `.log` run log are always written. The rest are opt-in:

| Flag | Adds |
| --- | --- |
| `--top` | the topology (`.top`) and the force-field include (`isolf.itp`) |
| `--inp` | GENESIS control files for a 3-step run (implies `--top`; also enabled by `--min`/`--equil`/`--prod`) |
| `--vmd` | a VMD script and `.psf` (also enabled by any `--vmd-*` flag) |
| `--gif` | a GIF animating the leaflet relaxation (also enabled by any `--gif-*` flag) |
| `--pdb` `--psf` `--crd` `--cif` | extra coordinate and structure formats |

The full file table and the reproducible run log are in [EXAMPLES.md → Output files](EXAMPLES.md#output-files).

## 🔬 Run in GENESIS

Run the control files in order. Each one continues from the previous step's restart:

```bash
atdyn min.inp     # 1. minimization
atdyn npt.inp     # 2. equilibration  (nvt.inp for a vesicle or tubule)
atdyn pro.inp     # 3. production
```

Run `isolf membrane --help`, `isolf vesicle --help`, or `isolf tubule --help` to see every option and its default.

## 🙏 Acknowledgements

Thanks to Hideto Tsubouchi for testing the code.

## 📄 License

Licensed under either of [Apache-2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT), at your option. Any contribution you submit is dual licensed as above, unless you say otherwise.

---

> Disclosure: This codebase was built with the help of a coding AI from my Python and Julia scripts, and is currently maintained by humans 🧑‍🔬 and ai 🤖.
