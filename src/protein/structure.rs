//! All-atom PDB and mmCIF parsers for coarse-graining protein structures.
//!
//! Coordinates are stored in Ångström (the native unit of PDB and mmCIF files)
//! so AICG2+ atomic contact classification matches `genesis_cg_tool` directly.

use std::fs;
use std::path::Path;

use crate::error::BuildError;
use crate::protein::params::canonical_aa;

/// One all-atom entry parsed from a PDB or mmCIF file.
#[derive(Debug, Clone)]
pub struct AllAtom {
    /// Trimmed PDB/mmCIF atom name (e.g. `"CA"`, `"NZ"`, `"OE1"`).
    pub name: String,
    /// Cartesian position in Ångström.
    pub pos: [f64; 3],
}

/// One protein residue with its constituent atoms (in input order).
#[derive(Debug, Clone)]
pub struct AllAtomResidue {
    /// Canonical 3-letter amino-acid name (e.g. `"HIS"` for `HSD`/`HSE`/`HSP`).
    pub name: String,
    /// Original residue sequence number from the input file.
    pub seq_id: i32,
    /// Indices into [`AllAtomStructure::atoms`].
    pub atoms: Vec<usize>,
}

/// One contiguous protein polypeptide chain.
#[derive(Debug, Clone)]
pub struct AllAtomChain {
    /// PDB chain identifier (e.g. `'A'`), or `' '` if absent.
    pub chain_id: char,
    /// Segment identifier (columns 73–76 in PDB, if present).
    pub seg_name: String,
    /// Indices into [`AllAtomStructure::residues`].
    pub residues: Vec<usize>,
}

/// All-atom protein structure ready for coarse-graining.
#[derive(Debug, Clone)]
pub struct AllAtomStructure {
    /// Every parsed atom across all chains.
    pub atoms: Vec<AllAtom>,
    /// Protein residues.
    pub residues: Vec<AllAtomResidue>,
    /// Protein chains.
    pub chains: Vec<AllAtomChain>,
}

/// Read an all-atom protein structure from a `.pdb`, `.ent`, `.cif`, or `.mmcif` file.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] on I/O or syntax errors, or
/// [`BuildError::EmptyProtein`] if no protein residues are found.
pub fn read_structure(path: &Path) -> Result<AllAtomStructure, BuildError> {
    let content = fs::read_to_string(path).map_err(|e| {
        BuildError::StructureParse(format!("failed to read {}: {e}", path.display()))
    })?;
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_ascii_lowercase());
    let is_cif = matches!(ext.as_deref(), Some("cif" | "mmcif"))
        || content.lines().any(|l| l.starts_with("_atom_site."));
    if is_cif {
        parse_cif(&content)
    } else {
        parse_pdb(&content)
    }
}

/// Slice fixed-width ASCII columns `[start..end]` safely (returning `""` if the
/// line is shorter).
fn col(line: &str, start: usize, end: usize) -> &str {
    let bytes = line.as_bytes();
    if start >= bytes.len() {
        return "";
    }
    let hi = end.min(bytes.len());
    std::str::from_utf8(&bytes[start..hi]).unwrap_or("")
}

/// Flush a pending residue into `residues` + `chain_residues` if it is a protein residue.
fn flush_residue(
    raw_res_name: &str,
    seq_id: i32,
    res_atoms: &mut Vec<usize>,
    residues: &mut Vec<AllAtomResidue>,
    chain_residues: &mut Vec<usize>,
) {
    if res_atoms.is_empty() {
        return;
    }
    let atoms = std::mem::take(res_atoms);
    if let Some(canon) = canonical_aa(raw_res_name) {
        let idx = residues.len();
        residues.push(AllAtomResidue {
            name: canon.to_string(),
            seq_id,
            atoms,
        });
        chain_residues.push(idx);
    }
}

/// Flush a pending chain into `chains` if it holds at least one protein residue.
fn flush_chain(
    chain_id: char,
    seg_name: &str,
    chain_residues: &mut Vec<usize>,
    chains: &mut Vec<AllAtomChain>,
) {
    if chain_residues.is_empty() {
        return;
    }
    let residues = std::mem::take(chain_residues);
    chains.push(AllAtomChain {
        chain_id,
        seg_name: seg_name.to_string(),
        residues,
    });
}

/// Parse an all-atom PDB string (first model only) into protein chains.
///
/// Follows PDB fixed-column conventions and `genesis_cg_tool`'s chain/residue
/// grouping: `TER` / `END` / `ENDMDL` terminate a chain, chain-ID changes also
/// start a new chain, alternate conformations other than `' '`, `'A'`, `'1'`
/// are skipped, and non-protein residues are ignored.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] on malformed coordinates or residue
/// numbers, or [`BuildError::EmptyProtein`] when no protein residues exist.
pub fn parse_pdb(content: &str) -> Result<AllAtomStructure, BuildError> {
    let mut atoms = Vec::new();
    let mut residues = Vec::new();
    let mut chains = Vec::new();

    let mut cur_chain_id: Option<char> = None;
    let mut cur_seg_name = String::new();
    let mut cur_res_key: Option<(i32, char)> = None;
    let mut cur_res_name = String::new();
    let mut cur_res_atoms: Vec<usize> = Vec::new();
    let mut cur_chain_res: Vec<usize> = Vec::new();

    for (line_no, line) in content.lines().enumerate() {
        if line.starts_with("ENDMDL") {
            flush_residue(
                &cur_res_name,
                cur_res_key.map_or(0, |k| k.0),
                &mut cur_res_atoms,
                &mut residues,
                &mut cur_chain_res,
            );
            flush_chain(
                cur_chain_id.unwrap_or(' '),
                &cur_seg_name,
                &mut cur_chain_res,
                &mut chains,
            );
            break;
        }
        if line.starts_with("TER") || line.starts_with("END") {
            flush_residue(
                &cur_res_name,
                cur_res_key.map_or(0, |k| k.0),
                &mut cur_res_atoms,
                &mut residues,
                &mut cur_chain_res,
            );
            flush_chain(
                cur_chain_id.unwrap_or(' '),
                &cur_seg_name,
                &mut cur_chain_res,
                &mut chains,
            );
            cur_chain_id = None;
            cur_res_key = None;
            continue;
        }
        if !line.starts_with("ATOM  ") && !line.starts_with("ATOM\t") {
            continue;
        }

        // Altloc (column 17, 0-based index 16): keep blank, 'A', or '1'.
        let alt = col(line, 16, 17).chars().next().unwrap_or(' ');
        if !matches!(alt, ' ' | 'A' | '1') {
            continue;
        }

        let atom_name = col(line, 12, 16).trim();
        let res_name = col(line, 17, 21).trim();
        let chain_id = col(line, 21, 22).chars().next().unwrap_or(' ');
        let res_seq_str = col(line, 22, 26).trim();
        let res_seq: i32 = res_seq_str.parse().map_err(|_| {
            BuildError::StructureParse(format!(
                "PDB line {}: invalid residue number '{res_seq_str}'",
                line_no + 1
            ))
        })?;
        let ins_code = col(line, 26, 27).chars().next().unwrap_or(' ');

        let parse_coord = |start: usize, end: usize, axis: &str| -> Result<f64, BuildError> {
            let s = col(line, start, end).trim();
            s.parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| {
                    BuildError::StructureParse(format!(
                        "PDB line {}: invalid {axis} coordinate '{s}'",
                        line_no + 1
                    ))
                })
        };
        let x = parse_coord(30, 38, "x")?;
        let y = parse_coord(38, 46, "y")?;
        let z = parse_coord(46, 54, "z")?;
        // Columns 67–76 (matching `genesis_cg_tool` line 64) hold segid in CHARMM/PDB.
        let seg_name = col(line, 72, 76).trim();

        // Start a new chain if the chain ID changed without an explicit TER.
        if let Some(prev_ch) = cur_chain_id
            && chain_id != prev_ch
        {
            flush_residue(
                &cur_res_name,
                cur_res_key.map_or(0, |k| k.0),
                &mut cur_res_atoms,
                &mut residues,
                &mut cur_chain_res,
            );
            flush_chain(prev_ch, &cur_seg_name, &mut cur_chain_res, &mut chains);
            cur_res_key = None;
        }

        let res_key = (res_seq, ins_code);
        if cur_res_key != Some(res_key) {
            flush_residue(
                &cur_res_name,
                cur_res_key.map_or(0, |k| k.0),
                &mut cur_res_atoms,
                &mut residues,
                &mut cur_chain_res,
            );
            cur_res_key = Some(res_key);
            cur_res_name = res_name.to_string();
        }
        cur_chain_id = Some(chain_id);
        cur_seg_name = seg_name.to_string();

        let atom_idx = atoms.len();
        atoms.push(AllAtom {
            name: atom_name.to_string(),
            pos: [x, y, z],
        });
        cur_res_atoms.push(atom_idx);
    }

    flush_residue(
        &cur_res_name,
        cur_res_key.map_or(0, |k| k.0),
        &mut cur_res_atoms,
        &mut residues,
        &mut cur_chain_res,
    );
    flush_chain(
        cur_chain_id.unwrap_or(' '),
        &cur_seg_name,
        &mut cur_chain_res,
        &mut chains,
    );

    if chains.is_empty() {
        return Err(BuildError::EmptyProtein);
    }
    Ok(AllAtomStructure {
        atoms,
        residues,
        chains,
    })
}

/// Tokenize one mmCIF data line, honoring `'...'` and `"..."` quoted strings.
fn tokenize_cif_line(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let n = bytes.len();
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        while i < n && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= n || bytes[i] == b'#' {
            break;
        }
        if bytes[i] == b'\'' || bytes[i] == b'"' {
            let q = bytes[i];
            i += 1;
            let start = i;
            while i < n && !(bytes[i] == q && (i + 1 == n || bytes[i + 1].is_ascii_whitespace())) {
                i += 1;
            }
            out.push(&line[start..i]);
            if i < n {
                i += 1;
            }
        } else {
            let start = i;
            while i < n && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push(&line[start..i]);
        }
    }
    out
}

/// Parse an all-atom mmCIF string (`_atom_site.*` loop, first model only) into
/// protein chains.
///
/// # Errors
///
/// Returns [`BuildError::StructureParse`] if required `_atom_site` columns are
/// missing or malformed, or [`BuildError::EmptyProtein`] if no protein residues
/// are present.
struct AtomSiteCols {
    group: Option<usize>,
    atom: usize,
    alt: Option<usize>,
    comp: usize,
    asym: Option<usize>,
    seq: usize,
    ins: Option<usize>,
    x: usize,
    y: usize,
    z: usize,
    model: Option<usize>,
}

impl AtomSiteCols {
    fn from_headers(headers: &[String]) -> Result<Self, BuildError> {
        let find = |k: &str| headers.iter().position(|h| h == k);
        let req = |opt: Option<usize>, label: &str| {
            opt.ok_or_else(|| {
                BuildError::StructureParse(format!(
                    "mmCIF _atom_site loop is missing required column {label}"
                ))
            })
        };
        Ok(Self {
            group: find("_atom_site.group_PDB"),
            atom: req(
                find("_atom_site.label_atom_id").or_else(|| find("_atom_site.auth_atom_id")),
                "label_atom_id/auth_atom_id",
            )?,
            alt: find("_atom_site.label_alt_id"),
            comp: req(
                find("_atom_site.label_comp_id").or_else(|| find("_atom_site.auth_comp_id")),
                "label_comp_id/auth_comp_id",
            )?,
            asym: find("_atom_site.label_asym_id").or_else(|| find("_atom_site.auth_asym_id")),
            seq: req(
                find("_atom_site.label_seq_id").or_else(|| find("_atom_site.auth_seq_id")),
                "label_seq_id/auth_seq_id",
            )?,
            ins: find("_atom_site.pdbx_PDB_ins_code"),
            x: req(find("_atom_site.Cartn_x"), "Cartn_x")?,
            y: req(find("_atom_site.Cartn_y"), "Cartn_y")?,
            z: req(find("_atom_site.Cartn_z"), "Cartn_z")?,
            model: find("_atom_site.pdbx_PDB_model_num"),
        })
    }
}

pub fn parse_cif(content: &str) -> Result<AllAtomStructure, BuildError> {
    let mut in_loop = false;
    let mut reading_headers = false;
    let mut headers: Vec<String> = Vec::new();
    let mut in_multiline = false;
    let mut cols: Option<AtomSiteCols> = None;

    let mut atoms = Vec::new();
    let mut residues = Vec::new();
    let mut chains = Vec::new();

    let mut first_model: Option<String> = None;
    let mut cur_asym: Option<String> = None;
    let mut cur_seq: Option<i32> = None;
    let mut cur_ins: Option<String> = None;
    let mut cur_res_name = String::new();
    let mut cur_res_atoms: Vec<usize> = Vec::new();
    let mut cur_chain_res: Vec<usize> = Vec::new();
    let mut found_atom_site = false;

    for (line_no, line) in content.lines().enumerate() {
        if line.starts_with(';') {
            in_multiline = !in_multiline;
            continue;
        }
        if in_multiline {
            continue;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            if cols.is_some() && !reading_headers {
                cols = None;
                in_loop = false;
            }
            continue;
        }
        if trimmed == "loop_" {
            in_loop = true;
            reading_headers = true;
            cols = None;
            headers.clear();
            continue;
        }
        if trimmed.starts_with("data_") || trimmed.starts_with("save_") || trimmed == "stop_" {
            in_loop = false;
            reading_headers = false;
            cols = None;
            continue;
        }

        if in_loop && reading_headers {
            if trimmed.starts_with('_') && !trimmed.contains(char::is_whitespace) {
                headers.push(trimmed.to_string());
                continue;
            }
            reading_headers = false;
            if headers.iter().any(|h| h.starts_with("_atom_site.")) {
                found_atom_site = true;
                cols = Some(AtomSiteCols::from_headers(&headers)?);
            }
        }

        let Some(c) = cols.as_ref() else {
            continue;
        };
        if trimmed.starts_with('_') {
            cols = None;
            in_loop = false;
            continue;
        }

        let tokens = tokenize_cif_line(trimmed);
        if tokens.len() < headers.len() {
            continue;
        }
        if let Some(g) = c.group
            && tokens[g] != "ATOM"
        {
            continue;
        }
        if let Some(m) = c.model {
            let mv = tokens[m];
            match &first_model {
                None => first_model = Some(mv.to_string()),
                Some(first) if first != mv => continue,
                _ => {}
            }
        }
        if let Some(a) = c.alt {
            let av = tokens[a];
            if !matches!(av, "." | "?" | "A" | "1") {
                continue;
            }
        }

        let atom_name = tokens[c.atom];
        let res_name = tokens[c.comp];
        let asym = c.asym.map_or("A", |i| tokens[i]);
        let seq_str = tokens[c.seq];
        let seq_id: i32 = match seq_str.parse() {
            Ok(v) => v,
            Err(_) => continue, // skip `.` or `?` non-polymer records
        };
        let ins = c.ins.map_or(".", |i| tokens[i]);

        let parse_f = |idx: usize, axis: &str| -> Result<f64, BuildError> {
            let tok = tokens[idx];
            tok.parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| {
                    BuildError::StructureParse(format!(
                        "mmCIF line {}: invalid {axis} coordinate '{tok}'",
                        line_no + 1
                    ))
                })
        };
        let x = parse_f(c.x, "x")?;
        let y = parse_f(c.y, "y")?;
        let z = parse_f(c.z, "z")?;

        let chain_break = cur_asym.as_deref().is_some_and(|prev| prev != asym)
            || cur_seq.is_some_and(|prev| seq_id > prev + 1);
        if chain_break {
            flush_residue(
                &cur_res_name,
                cur_seq.unwrap_or(0),
                &mut cur_res_atoms,
                &mut residues,
                &mut cur_chain_res,
            );
            let ch = cur_asym
                .as_deref()
                .and_then(|s| s.chars().next())
                .unwrap_or(' ');
            flush_chain(ch, "", &mut cur_chain_res, &mut chains);
            cur_seq = None;
            cur_ins = None;
        }

        if cur_seq != Some(seq_id) || cur_ins.as_deref() != Some(ins) {
            flush_residue(
                &cur_res_name,
                cur_seq.unwrap_or(0),
                &mut cur_res_atoms,
                &mut residues,
                &mut cur_chain_res,
            );
            cur_seq = Some(seq_id);
            cur_ins = Some(ins.to_string());
            cur_res_name = res_name.to_string();
        }
        cur_asym = Some(asym.to_string());

        let atom_idx = atoms.len();
        atoms.push(AllAtom {
            name: atom_name.to_string(),
            pos: [x, y, z],
        });
        cur_res_atoms.push(atom_idx);
    }

    flush_residue(
        &cur_res_name,
        cur_seq.unwrap_or(0),
        &mut cur_res_atoms,
        &mut residues,
        &mut cur_chain_res,
    );
    let ch = cur_asym
        .as_deref()
        .and_then(|s| s.chars().next())
        .unwrap_or(' ');
    flush_chain(ch, "", &mut cur_chain_res, &mut chains);

    if !found_atom_site {
        return Err(BuildError::StructureParse(
            "mmCIF file contains no _atom_site loop".to_string(),
        ));
    }
    if chains.is_empty() {
        return Err(BuildError::EmptyProtein);
    }
    Ok(AllAtomStructure {
        atoms,
        residues,
        chains,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_PDB: &str = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  N   HSD A   2       2.000   1.200   0.000  1.00  0.00           N
ATOM      4  CA  HSD A   2       3.400   1.200   0.000  1.00  0.00           C
TER
ATOM      5  N   GLY B   1      10.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  GLY B   1      11.400   0.000   0.000  1.00  0.00           C
HETATM    7  O   HOH C   1      20.000   0.000   0.000  1.00  0.00           O
END
";

    const SAMPLE_CIF: &str = "\
data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N . ALA A 1 0.0 0.0 0.0 1
ATOM 2 CA A ALA A 1 1.458 0.0 0.0 1
ATOM 3 CA B ALA A 1 1.999 0.0 0.0 1
ATOM 4 N . HSD A 2 2.0 1.2 0.0 1
ATOM 5 CA . HSD A 2 3.4 1.2 0.0 1
ATOM 6 N . GLY B 1 10.0 0.0 0.0 1
ATOM 7 CA . GLY B 1 11.4 0.0 0.0 1
ATOM 8 CA . GLY B 1 99.0 0.0 0.0 2
#
";

    #[test]
    fn parses_pdb_chains_and_aliases() {
        let s = parse_pdb(SAMPLE_PDB).unwrap();
        assert_eq!(s.chains.len(), 2);
        assert_eq!(s.chains[0].chain_id, 'A');
        assert_eq!(s.chains[1].chain_id, 'B');
        assert_eq!(s.residues.len(), 3);
        assert_eq!(s.residues[0].name, "ALA");
        assert_eq!(s.residues[1].name, "HIS"); // HSD -> HIS
        assert_eq!(s.residues[2].name, "GLY");
    }

    #[test]
    fn parses_cif_filtering_altloc_and_model_2() {
        let s = parse_cif(SAMPLE_CIF).unwrap();
        assert_eq!(s.chains.len(), 2);
        assert_eq!(s.residues.len(), 3);
        assert_eq!(s.atoms.len(), 6); // altloc B and model 2 skipped
        assert_eq!(s.residues[1].name, "HIS");
    }

    #[test]
    fn empty_pdb_returns_empty_protein() {
        assert_eq!(parse_pdb("END\n").unwrap_err(), BuildError::EmptyProtein);
    }
}
