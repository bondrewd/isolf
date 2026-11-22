//! Verbatim AICG2+ and Clementi Gō force-field parameters for GENESIS RESIDCG.
//!
//! All values match `genesis_cg_tool` (`param/atom_types.itp`,
//! `param/flexible_local_angle.itp`, `param/flexible_local_dihedral.itp`, and
//! `src/lib/constants.jl`). Energies in the rendered `aicg2p.itp` are in
//! kJ/mol and lengths are in nm.

use std::fmt::Write;

/// kcal/mol → kJ/mol conversion factor.
pub const CAL_TO_KJ: f64 = 4.184;

/// AICG2+ covalent Cα–Cα bond force constant (kcal mol⁻¹ Å⁻²).
pub const AICG_BOND_K: f64 = 110.40;
/// AICG2+ Gaussian 1–3 distance width (Å).
pub const AICG_13_SIGMA: f64 = 0.15;
/// AICG2+ Gaussian 1–4 dihedral width (rad).
pub const AICG_14_SIGMA: f64 = 0.15;
/// Atomic heavy-atom contact cutoff for Gō native contacts (Å).
pub const AICG_GO_ATOMIC_CUTOFF: f64 = 6.5;
/// Short-range atomic contact cutoff for AICG2+ energy classification (Å).
pub const AICG_ATOMIC_CUTOFF: f64 = 5.0;
/// Hydrogen-bond atomic cutoff (Å).
pub const AICG_HYDROGEN_BOND_CUTOFF: f64 = 3.2;
/// Salt-bridge atomic cutoff (Å).
pub const AICG_SALT_BRIDGE_CUTOFF: f64 = 3.5;
/// Upper clamp for raw AICG2+ local/contact energy (kcal/mol).
pub const AICG_ENE_UPPER_LIM: f64 = -0.5;
/// Lower clamp for raw AICG2+ local/contact energy (kcal/mol).
pub const AICG_ENE_LOWER_LIM: f64 = -5.0;

/// Average scaling targets (kcal/mol, `aicg-scale = 0`).
pub const AICG_13_AVE: f64 = 1.72;
pub const AICG_14_AVE: f64 = 1.23;
pub const AICG_CONTACT_AVE: f64 = 0.55;

/// General scaling multipliers (`aicg-scale = 1`, default).
pub const AICG_13_GEN: f64 = 1.11;
pub const AICG_14_GEN: f64 = 0.87;
pub const AICG_CONTACT_GEN: f64 = 0.32;

/// Clementi Gō energy unit and force constants (kcal/mol, Å, rad).
pub const CCGO_EPSILON: f64 = 1.0;
pub const CCGO_BOND_K: f64 = 100.00 * CCGO_EPSILON;
pub const CCGO_ANGL_K: f64 = 20.00 * CCGO_EPSILON;
pub const CCGO_DIHE_K_1: f64 = CCGO_EPSILON;
pub const CCGO_DIHE_K_3: f64 = 0.5 * CCGO_EPSILON;
pub const CCGO_NATIVE_EPSILON: f64 = CCGO_EPSILON;

/// 17-class AICG2+ atomic pairwise energy weights (kcal/mol), ordered to match
/// `genesis_cg_tool` (`BB_HB..OFFST`).
pub const AICG_PAIRWISE_ENERGY: [f64; 17] = [
    -1.4247,  // 0:  BB_HB (backbone–backbone H-bond)
    -0.4921,  // 1:  BB_DA (backbone–backbone donor–acceptor)
    -0.2404,  // 2:  BB_CX (backbone–backbone carbon–X)
    -0.1035,  // 3:  BB_XX (backbone–backbone other)
    -5.7267,  // 4:  SS_HB (sidechain–sidechain H-bond)
    -12.4878, // 5:  SS_SB (sidechain–sidechain salt bridge)
    -0.0308,  // 6:  SS_DA (sidechain–sidechain donor–acceptor)
    -0.1113,  // 7:  SS_CX (sidechain–sidechain carbon–X)
    -0.2168,  // 8:  SS_QX (sidechain–sidechain charge–X)
    0.2306,   // 9:  SS_XX (sidechain–sidechain other)
    -3.4819,  // 10: SB_HB (sidechain–backbone H-bond)
    -0.1809,  // 11: SB_DA (sidechain–backbone donor–acceptor)
    -0.1209,  // 12: SB_CX (sidechain–backbone carbon–X)
    -0.2984,  // 13: SB_QX (sidechain–backbone charge–X)
    -0.0487,  // 14: SB_XX (sidechain–backbone other)
    -0.0395,  // 15: LR_CT (long-range 5.0–6.5 Å contacts)
    -0.1051,  // 16: OFFST (constant offset)
];

/// Indices into [`AICG_PAIRWISE_ENERGY`].
pub const ITYPE_BB_HB: usize = 0;
pub const ITYPE_BB_DA: usize = 1;
pub const ITYPE_BB_CX: usize = 2;
pub const ITYPE_BB_XX: usize = 3;
pub const ITYPE_SS_HB: usize = 4;
pub const ITYPE_SS_SB: usize = 5;
pub const ITYPE_SS_DA: usize = 6;
pub const ITYPE_SS_CX: usize = 7;
pub const ITYPE_SS_QX: usize = 8;
pub const ITYPE_SS_XX: usize = 9;
pub const ITYPE_SB_HB: usize = 10;
pub const ITYPE_SB_DA: usize = 11;
pub const ITYPE_SB_CX: usize = 12;
pub const ITYPE_SB_QX: usize = 13;
pub const ITYPE_SB_XX: usize = 14;
pub const ITYPE_LR_CT: usize = 15;
pub const ITYPE_OFFST: usize = 16;

/// Static AICG2+ residue bead parameters (from `atom_types.itp`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AminoAcid {
    /// Standard 3-letter residue code (e.g. `"ALA"`).
    pub name: &'static str,
    /// 1-letter FASTA code.
    pub short: char,
    /// Atom type `n` column in `atom_types.itp`.
    pub n_val: u8,
    /// Residue mass (g/mol).
    pub mass: f64,
    /// Standard net charge (e).
    pub charge: f64,
    /// Excluded-volume diameter `rmin` (nm).
    pub rmin: f64,
    /// Excluded-volume epsilon `eps` (kJ/mol).
    pub eps: f64,
}

/// The 20 standard amino acids in canonical alphabetical order.
#[allow(clippy::approx_constant)]
pub const AMINO_ACIDS: [AminoAcid; 20] = [
    AminoAcid {
        name: "ALA",
        short: 'A',
        n_val: 6,
        mass: 71.090,
        charge: 0.0,
        rmin: 0.3943,
        eps: 2.5104,
    },
    AminoAcid {
        name: "ARG",
        short: 'R',
        n_val: 6,
        mass: 156.190,
        charge: 1.0,
        rmin: 0.4166,
        eps: 2.5104,
    },
    AminoAcid {
        name: "ASN",
        short: 'N',
        n_val: 6,
        mass: 114.110,
        charge: 0.0,
        rmin: 0.4328,
        eps: 2.5104,
    },
    AminoAcid {
        name: "ASP",
        short: 'D',
        n_val: 6,
        mass: 115.090,
        charge: -1.0,
        rmin: 0.4419,
        eps: 2.5104,
    },
    AminoAcid {
        name: "CYS",
        short: 'C',
        n_val: 6,
        mass: 103.150,
        charge: 0.0,
        rmin: 0.4449,
        eps: 2.5104,
    },
    AminoAcid {
        name: "GLN",
        short: 'Q',
        n_val: 6,
        mass: 128.140,
        charge: 0.0,
        rmin: 0.3955,
        eps: 2.5104,
    },
    AminoAcid {
        name: "GLU",
        short: 'E',
        n_val: 6,
        mass: 129.120,
        charge: -1.0,
        rmin: 0.4316,
        eps: 2.5104,
    },
    AminoAcid {
        name: "GLY",
        short: 'G',
        n_val: 6,
        mass: 57.050,
        charge: 0.0,
        rmin: 0.3448,
        eps: 2.5104,
    },
    AminoAcid {
        name: "HIS",
        short: 'H',
        n_val: 6,
        mass: 137.140,
        charge: 0.0,
        rmin: 0.4696,
        eps: 2.5104,
    },
    AminoAcid {
        name: "ILE",
        short: 'I',
        n_val: 6,
        mass: 113.160,
        charge: 0.0,
        rmin: 0.4352,
        eps: 2.5104,
    },
    AminoAcid {
        name: "LEU",
        short: 'L',
        n_val: 6,
        mass: 113.160,
        charge: 0.0,
        rmin: 0.4342,
        eps: 2.5104,
    },
    AminoAcid {
        name: "LYS",
        short: 'K',
        n_val: 6,
        mass: 128.170,
        charge: 1.0,
        rmin: 0.4530,
        eps: 2.5104,
    },
    AminoAcid {
        name: "MET",
        short: 'M',
        n_val: 7,
        mass: 131.190,
        charge: 0.0,
        rmin: 0.4424,
        eps: 2.5104,
    },
    AminoAcid {
        name: "PHE",
        short: 'F',
        n_val: 6,
        mass: 147.180,
        charge: 0.0,
        rmin: 0.4403,
        eps: 2.5104,
    },
    AminoAcid {
        name: "PRO",
        short: 'P',
        n_val: 6,
        mass: 97.120,
        charge: 0.0,
        rmin: 0.4046,
        eps: 2.5104,
    },
    AminoAcid {
        name: "SER",
        short: 'S',
        n_val: 6,
        mass: 87.080,
        charge: 0.0,
        rmin: 0.3927,
        eps: 2.5104,
    },
    AminoAcid {
        name: "THR",
        short: 'T',
        n_val: 6,
        mass: 101.110,
        charge: 0.0,
        rmin: 0.4036,
        eps: 2.5104,
    },
    AminoAcid {
        name: "TRP",
        short: 'W',
        n_val: 6,
        mass: 186.210,
        charge: 0.0,
        rmin: 0.4476,
        eps: 2.5104,
    },
    AminoAcid {
        name: "TYR",
        short: 'Y',
        n_val: 6,
        mass: 163.180,
        charge: 0.0,
        rmin: 0.4135,
        eps: 2.5104,
    },
    AminoAcid {
        name: "VAL",
        short: 'V',
        n_val: 6,
        mass: 99.140,
        charge: 0.0,
        rmin: 0.4440,
        eps: 2.5104,
    },
];

/// Map a PDB/mmCIF protein residue name (including protonation aliases) to its
/// canonical 3-letter AICG2+ residue name, or `None` if not a protein residue.
pub fn canonical_aa(name: &str) -> Option<&'static str> {
    match name.trim() {
        "ALA" => Some("ALA"),
        "ARG" => Some("ARG"),
        "ASN" => Some("ASN"),
        "ASP" => Some("ASP"),
        "CYS" | "CYM" => Some("CYS"),
        "GLN" => Some("GLN"),
        "GLU" => Some("GLU"),
        "GLY" => Some("GLY"),
        "HIS" | "HSD" | "HSE" | "HSP" | "HID" | "HIE" | "HIP" => Some("HIS"),
        "ILE" => Some("ILE"),
        "LEU" => Some("LEU"),
        "LYS" => Some("LYS"),
        "MET" => Some("MET"),
        "PHE" => Some("PHE"),
        "PRO" => Some("PRO"),
        "SER" => Some("SER"),
        "THR" => Some("THR"),
        "TRP" => Some("TRP"),
        "TYR" => Some("TYR"),
        "VAL" => Some("VAL"),
        "UNK" => Some("UNK"),
        _ => None,
    }
}

/// Look up AICG2+ bead parameters for a canonical (or alias) amino acid name.
pub fn aa_params(name: &str) -> Option<&'static AminoAcid> {
    let canon = canonical_aa(name)?;
    AMINO_ACIDS.iter().find(|aa| aa.name == canon)
}

/// Mass (g/mol) for a canonical amino acid (or `UNK` = 100.0).
pub fn aa_mass(canonical: &str) -> f64 {
    aa_params(canonical).map_or(100.0, |aa| aa.mass)
}

/// Standard charge (e) for a canonical amino acid (or `UNK` = 0.0).
pub fn aa_charge(canonical: &str) -> f64 {
    aa_params(canonical).map_or(0.0, |aa| aa.charge)
}

/// Excluded-volume diameter `rmin` (nm) for a canonical amino acid (or `UNK` = 0.42).
pub fn aa_rmin(canonical: &str) -> f64 {
    aa_params(canonical).map_or(0.42, |aa| aa.rmin)
}

/// Shared 10 spline-angle knots (rad) for `[ flexible_local_angle ]`.
const FLEX_ANGLE_X: [f64; 10] = [
    1.3090, 1.4835, 1.6581, 1.8326, 2.0071, 2.1817, 2.3562, 2.5307, 2.7053, 2.8798,
];

/// Verbatim `(y, y2)` cubic-spline values for `[ flexible_local_angle ]`,
/// indexed in the same order as [`AMINO_ACIDS`].
const FLEX_ANGLE_Y_Y2: [[(f64, f64); 10]; 20] = [
    // ALA
    [
        (2.0920E+01, 0.0000E+00),
        (5.6066E+00, 6.3580E+02),
        (3.5146E+00, 6.1128E+01),
        (4.8953E+00, -1.9619E+02),
        (3.4309E+00, 1.6334E+02),
        (4.1840E+00, -2.0334E+01),
        (5.3137E+00, -7.7822E+00),
        (6.3597E+00, 3.5062E+01),
        (1.3389E+01, 1.0461E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // ARG
    [
        (2.0920E+01, 0.0000E+00),
        (5.9831E+00, 5.9948E+02),
        (3.5146E+00, 5.8074E+01),
        (4.5606E+00, -1.3949E+02),
        (3.6819E+00, 1.2083E+02),
        (4.2258E+00, -6.3639E+01),
        (4.6442E+00, 1.0904E+02),
        (6.6944E+00, -5.1045E+01),
        (1.3640E+01, 1.0594E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // ASN
    [
        (2.0920E+01, 0.0000E+00),
        (5.7739E+00, 5.7296E+02),
        (2.6778E+00, 8.1797E+01),
        (3.8493E+00, -5.9455E+01),
        (4.4350E+00, 4.0668E+01),
        (5.3137E+00, -4.5564E+01),
        (5.8158E+00, 6.7404E+01),
        (7.1546E+00, -5.9162E+01),
        (1.3305E+01, 1.1171E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // ASP
    [
        (2.0920E+01, 0.0000E+00),
        (5.5229E+00, 6.0015E+02),
        (2.6778E+00, 7.1756E+01),
        (3.9330E+00, -7.9580E+01),
        (4.1840E+00, 4.8702E+01),
        (5.1882E+00, 3.3095E+01),
        (6.7362E+00, -7.3931E+01),
        (7.5730E+00, 1.2247E+02),
        (1.5146E+01, 9.1098E+02),
        (4.1840E+01, 0.0000E+00),
    ],
    // CYS
    [
        (2.0920E+01, 0.0000E+00),
        (6.6526E+00, 5.6957E+02),
        (4.0166E+00, 1.2761E+01),
        (4.1840E+00, -6.8534E+01),
        (3.4727E+00, 8.8241E+01),
        (4.0585E+00, -2.8870E+01),
        (4.6442E+00, 2.7322E+01),
        (6.0250E+00, 7.6232E+01),
        (1.3975E+01, 9.6169E+02),
        (4.1840E+01, 0.0000E+00),
    ],
    // GLN
    [
        (2.0920E+01, 0.0000E+00),
        (5.8994E+00, 5.9756E+02),
        (3.3054E+00, 5.7404E+01),
        (4.3932E+00, -1.0192E+02),
        (4.0166E+00, 6.1798E+01),
        (4.2258E+00, -2.9916E+01),
        (4.7279E+00, 1.1560E+02),
        (6.9454E+00, -9.4558E+01),
        (1.3431E+01, 1.1033E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // GLU
    [
        (2.0920E+01, 0.0000E+00),
        (5.4392E+00, 6.4044E+02),
        (3.2217E+00, 5.0836E+01),
        (4.5606E+00, -1.4322E+02),
        (3.8074E+00, 1.1000E+02),
        (4.3932E+00, -3.2970E+01),
        (5.2718E+00, 7.9622E+01),
        (7.3638E+00, -4.6484E+01),
        (1.4184E+01, 1.0377E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // GLY
    [
        (2.0920E+01, 0.0000E+00),
        (6.6944E+00, 4.1555E+02),
        (2.3430E+00, 2.8284E+02),
        (5.0626E+00, -1.5410E+02),
        (5.7739E+00, -6.1965E+01),
        (4.6861E+00, 4.7530E+01),
        (4.5187E+00, 5.3220E+01),
        (5.1463E+00, -1.0372E+02),
        (1.0544E+01, 1.3013E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // HIS
    [
        (2.0920E+01, 0.0000E+00),
        (6.2760E+00, 5.7262E+02),
        (3.4727E+00, 4.1840E+01),
        (4.1422E+00, -5.6024E+01),
        (4.1003E+00, 4.2091E+01),
        (4.3932E+00, -4.6401E+01),
        (4.4769E+00, 1.0234E+02),
        (5.9831E+00, -8.2760E+01),
        (1.2259E+01, 1.1683E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // ILE
    [
        (2.0920E+01, 0.0000E+00),
        (7.4894E+00, 5.4739E+02),
        (4.8953E+00, -5.5061E+01),
        (3.8493E+00, -2.2217E+01),
        (2.7196E+00, 1.2744E+02),
        (3.7238E+00, -6.7237E+01),
        (4.9371E+00, 1.8267E+02),
        (8.6190E+00, -1.7728E+02),
        (1.4979E+01, 1.0539E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // LEU
    [
        (2.0920E+01, 0.0000E+00),
        (6.4434E+00, 5.5605E+02),
        (3.6401E+00, 7.5145E+01),
        (4.4769E+00, -1.3966E+02),
        (3.3890E+00, 1.0443E+02),
        (3.6819E+00, -6.0668E+00),
        (5.1463E+00, 1.5054E+02),
        (8.8282E+00, -1.5937E+02),
        (1.5272E+01, 1.0309E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // LYS
    [
        (2.0920E+01, 0.0000E+00),
        (5.6902E+00, 6.2639E+02),
        (3.3472E+00, 3.2844E+01),
        (4.3514E+00, -9.8408E+01),
        (3.8911E+00, 7.2425E+01),
        (4.2677E+00, -2.6401E+01),
        (4.9371E+00, 9.0960E+01),
        (7.1965E+00, -2.4267E+01),
        (1.4477E+01, 9.9504E+02),
        (4.1840E+01, 0.0000E+00),
    ],
    // MET
    [
        (2.0920E+01, 0.0000E+00),
        (6.6107E+00, 5.3886E+02),
        (3.6401E+00, 7.7990E+01),
        (4.4350E+00, -1.0916E+02),
        (3.7238E+00, 6.2049E+01),
        (3.8074E+00, 1.7656E+01),
        (4.5606E+00, -7.5312E-01),
        (6.6944E+00, 2.5732E+02),
        (1.7238E+01, 6.2793E+02),
        (4.1840E+01, 0.0000E+00),
    ],
    // PHE
    [
        (2.0920E+01, 0.0000E+00),
        (6.9454E+00, 5.2313E+02),
        (3.8493E+00, 5.0208E+01),
        (4.1422E+00, -5.6358E+01),
        (3.8074E+00, 5.1631E+01),
        (4.0585E+00, -3.4769E+01),
        (4.2258E+00, 7.0961E+01),
        (5.8158E+00, 3.1087E+01),
        (1.3514E+01, 1.0080E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // PRO
    [
        (2.0920E+01, 0.0000E+00),
        (5.1882E+00, 6.6358E+02),
        (3.4309E+00, 9.8282E+01),
        (5.4392E+00, -3.1493E+02),
        (2.6778E+00, 2.2200E+02),
        (3.7238E+00, 1.7682E+02),
        (9.3303E+00, -3.1045E+01),
        (1.4351E+01, -1.6807E+02),
        (2.0041E+01, 8.3525E+02),
        (4.1840E+01, 0.0000E+00),
    ],
    // SER
    [
        (2.0920E+01, 0.0000E+00),
        (5.8158E+00, 6.0999E+02),
        (3.5146E+00, 8.1881E+01),
        (5.0626E+00, -1.7928E+02),
        (3.9748E+00, 1.1602E+02),
        (4.1422E+00, -3.7614E+01),
        (4.3932E+00, 5.0836E+01),
        (5.3974E+00, -1.7405E+01),
        (1.2092E+01, 1.1396E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // THR
    [
        (2.0920E+01, 0.0000E+00),
        (6.6526E+00, 5.2936E+02),
        (3.5146E+00, 7.4810E+01),
        (4.2258E+00, -7.0417E+01),
        (4.1422E+00, 5.0250E+01),
        (4.3514E+00, -7.2843E+01),
        (4.0166E+00, 1.3401E+02),
        (6.0250E+00, -1.6736E+00),
        (1.3765E+01, 1.0017E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // TRP
    [
        (2.0920E+01, 0.0000E+00),
        (6.1086E+00, 5.8911E+02),
        (3.5146E+00, 5.0041E+01),
        (4.3514E+00, -1.1351E+02),
        (3.6819E+00, 1.0740E+02),
        (4.3095E+00, -6.0584E+01),
        (4.7698E+00, 1.0201E+02),
        (6.7362E+00, -5.0836E+01),
        (1.3598E+01, 1.0655E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // TYR
    [
        (2.0920E+01, 0.0000E+00),
        (7.1546E+00, 4.9133E+02),
        (3.8074E+00, 8.6902E+01),
        (4.2677E+00, -8.8910E+01),
        (3.7656E+00, 7.9203E+01),
        (4.1003E+00, -6.3095E+01),
        (4.1003E+00, 1.0724E+02),
        (5.9831E+00, 4.9371E+00),
        (1.3640E+01, 1.0104E+03),
        (4.1840E+01, 0.0000E+00),
    ],
    // VAL
    [
        (2.0920E+01, 0.0000E+00),
        (7.3220E+00, 5.7338E+02),
        (4.9790E+00, -7.6484E+01),
        (3.9748E+00, -3.6401E+00),
        (3.0125E+00, 9.9203E+01),
        (3.6401E+00, -8.0082E+01),
        (4.2258E+00, 2.1288E+02),
        (7.9914E+00, -1.4523E+02),
        (1.5021E+01, 1.0109E+03),
        (4.1840E+01, 0.0000E+00),
    ],
];

/// Verbatim 7-term Fourier coefficients for `[ flexible_local_dihedral_angle ]`
/// when `AA2` is neither `GLY` nor `PRO`, indexed by `AA1` in [`AMINO_ACIDS`].
const FLEX_DIHEDRAL_BY_AA1: [[f64; 7]; 20] = [
    [
        9.2282E+00,
        9.1337E-01,
        -3.3263E-01,
        1.8870E-01,
        -1.3259E+00,
        6.9036E-02,
        -5.7530E-01,
    ], // ALA
    [
        9.2207E+00,
        9.8993E-01,
        -4.7823E-01,
        3.7865E-01,
        -1.2318E+00,
        -1.3180E-01,
        -5.3513E-01,
    ], // ARG
    [
        9.2546E+00,
        1.5272E+00,
        -7.8910E-01,
        5.0250E-01,
        -4.4267E-01,
        4.6442E-02,
        -4.5187E-01,
    ], // ASN
    [
        9.2031E+00,
        1.2117E+00,
        -7.5730E-01,
        2.8577E-01,
        -8.7655E-01,
        1.4309E-01,
        -3.6484E-01,
    ], // ASP
    [
        9.4634E+00,
        2.0167E+00,
        -9.9998E-02,
        2.2050E-01,
        -1.2104E+00,
        1.5146E-01,
        -1.7949E-01,
    ], // CYS
    [
        9.2307E+00,
        1.0878E+00,
        -5.4978E-01,
        4.7070E-01,
        -1.2401E+00,
        -2.1840E-01,
        -5.7153E-01,
    ], // GLN
    [
        9.2286E+00,
        6.1923E-01,
        -5.7321E-01,
        3.1631E-01,
        -1.4213E+00,
        -5.8994E-02,
        -5.3388E-01,
    ], // GLU
    [
        9.0860E+00,
        1.2970E-02,
        -8.7446E-02,
        8.7655E-01,
        3.1422E-01,
        -2.9372E-01,
        -3.6819E-01,
    ], // GLY
    [
        9.2157E+00,
        1.3560E+00,
        -4.8953E-01,
        2.6819E-01,
        -9.6609E-01,
        -1.2175E-01,
        -1.2259E-01,
    ], // HIS
    [
        9.5157E+00,
        1.9046E+00,
        -3.9622E-01,
        -4.9204E-01,
        -1.4481E+00,
        -9.5814E-02,
        -2.4602E-01,
    ], // ILE
    [
        9.3638E+00,
        1.5832E+00,
        -7.7781E-01,
        4.5606E-02,
        -1.3715E+00,
        -2.8284E-01,
        -3.2342E-01,
    ], // LEU
    [
        9.2115E+00,
        7.7948E-01,
        -6.5270E-01,
        3.8744E-01,
        -1.2853E+00,
        -1.9497E-01,
        -5.3011E-01,
    ], // LYS
    [
        9.2914E+00,
        1.4652E+00,
        -4.6777E-01,
        3.2719E-01,
        -1.1293E+00,
        -1.7071E-01,
        -3.2886E-01,
    ], // MET
    [
        9.3538E+00,
        1.7234E+00,
        -5.3681E-01,
        9.5814E-02,
        -1.1330E+00,
        -3.5564E-02,
        -2.6819E-01,
    ], // PHE
    [
        9.7525E+00,
        7.3053E-01,
        1.4393E-01,
        5.1463E-01,
        -2.6363E+00,
        3.7991E-01,
        -6.9622E-01,
    ], // PRO
    [
        9.2663E+00,
        1.0581E+00,
        -3.4644E-01,
        5.6442E-01,
        -1.4150E+00,
        7.6986E-02,
        -2.1840E-01,
    ], // SER
    [
        9.3065E+00,
        9.9161E-01,
        -3.4058E-01,
        8.0584E-01,
        -1.5104E+00,
        -2.6317E-01,
        -1.4979E-01,
    ], // THR
    [
        9.2889E+00,
        1.2661E+00,
        -7.7864E-01,
        -1.0878E-01,
        -1.2502E+00,
        4.2258E-02,
        -1.5481E-02,
    ], // TRP
    [
        9.3123E+00,
        1.6330E+00,
        -5.2174E-01,
        5.4392E-02,
        -1.0970E+00,
        1.2175E-01,
        -2.5815E-01,
    ], // TYR
    [
        9.4487E+00,
        1.7171E+00,
        -1.1213E-01,
        -2.2343E-01,
        -1.5217E+00,
        -5.8158E-02,
        -3.1464E-01,
    ], // VAL
];

/// Dihedral coefficients for `(*, GLY)`.
const FLEX_DIHEDRAL_ANY_GLY: [f64; 7] = [
    9.4006E+00,
    -6.8869E-01,
    4.6944E-01,
    1.5857E-01,
    -1.4841E+00,
    -6.9454E-01,
    -5.8367E-01,
];

/// Dihedral coefficients for `(GLY, PRO)`.
const FLEX_DIHEDRAL_GLY_PRO: [f64; 7] = [
    1.0417E+01,
    -1.8882E+00,
    3.5531E+00,
    -1.6736E-02,
    -4.6652E-01,
    1.0887E+00,
    2.6610E-01,
];

/// Dihedral coefficients for `(X != GLY, PRO)`.
const FLEX_DIHEDRAL_OTHER_PRO: [f64; 7] = [
    1.1338E+01,
    2.0292E+00,
    3.2217E+00,
    3.9664E-01,
    -3.1271E+00,
    7.8450E-01,
    -1.5477E+00,
];

/// Index of `GLY` and `PRO` in [`AMINO_ACIDS`].
const IDX_GLY: usize = 7;
const IDX_PRO: usize = 14;

/// Fixed-capacity stack buffer for zero-allocation scientific formatting.
struct StackBuf {
    buf: [u8; 32],
    len: usize,
}

impl std::fmt::Write for StackBuf {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        let bytes = s.as_bytes();
        let end = self.len + bytes.len();
        if end > self.buf.len() {
            return Err(std::fmt::Error);
        }
        self.buf[self.len..end].copy_from_slice(bytes);
        self.len = end;
        Ok(())
    }
}

/// Append a floating-point number in C/Fortran-style scientific notation with a
/// signed at-least-two-digit exponent (e.g. `3.8017E-01`, `-1.2656E+02`),
/// right-aligned to `width`, with zero heap allocations.
pub(crate) fn write_sci(out: &mut String, value: f64, width: usize) {
    let v = if value == 0.0 { 0.0 } else { value };
    let mut tmp = StackBuf {
        buf: [0; 32],
        len: 0,
    };
    let _ = write!(tmp, "{v:.4E}");
    let raw = std::str::from_utf8(&tmp.buf[..tmp.len]).unwrap_or("0.0000E0");
    let (mantissa, exp_str) = raw.split_once('E').unwrap_or((raw, "0"));
    let exp: i32 = exp_str.parse().unwrap_or(0);
    let abs_exp = exp.unsigned_abs();
    let exp_digits = if abs_exp >= 100 { 3 } else { 2 };
    let formatted_len = mantissa.len() + 2 + exp_digits;
    for _ in formatted_len..width {
        out.push(' ');
    }
    out.push_str(mantissa);
    out.push('E');
    out.push(if exp < 0 { '-' } else { '+' });
    let _ = write!(out, "{abs_exp:02}");
}

/// Render a self-contained GENESIS `aicg2p.itp` parameter file containing
/// `[ defaults ]`, protein `[ atomtypes ]` (20 AA + 20 `iAA`),
/// `[ flexible_local_angle ]`, and `[ flexible_local_dihedral_angle ]`.
pub fn render_aicg2p_itp() -> String {
    let mut s = String::with_capacity(96 * 1024);

    s.push_str("[ defaults ]\n");
    s.push_str(";nbfunc comb-rule gen-pairs fudge_lj fudge_qq\n");
    s.push_str("      1         2        no      1.0      1.0\n\n");

    s.push_str("[ atomtypes ]\n");
    s.push_str(";PROTEIN\n");
    s.push_str(";name  n    mass  charge ptype       rmin      eps\n");
    for (pad, prefix) in [("  ", ""), (" ", "i")] {
        for aa in &AMINO_ACIDS {
            let _ = writeln!(
                s,
                "{pad}{prefix}{} {:>2} {:>7.3} {:>7.3}     A {:>10.4} {:>8.4}",
                aa.name, aa.n_val, aa.mass, aa.charge, aa.rmin, aa.eps
            );
        }
    }
    s.push('\n');

    s.push_str("[ flexible_local_angle ]\n");
    for prefix in ["", "i"] {
        for (idx, aa) in AMINO_ACIDS.iter().enumerate() {
            let _ = writeln!(s, "{prefix}{} 1 10 \\", aa.name);
            for (k, (&x, &(y, y2))) in FLEX_ANGLE_X.iter().zip(&FLEX_ANGLE_Y_Y2[idx]).enumerate() {
                write_sci(&mut s, x, 15);
                write_sci(&mut s, y, 15);
                write_sci(&mut s, y2, 15);
                if k + 1 < 10 {
                    s.push_str(" \\\n");
                } else {
                    s.push('\n');
                }
            }
            s.push('\n');
        }
    }

    s.push_str("[ flexible_local_dihedral_angle ]\n");
    for prefix in ["", "i"] {
        for (i, aa1) in AMINO_ACIDS.iter().enumerate() {
            for (j, aa2) in AMINO_ACIDS.iter().enumerate() {
                let coeffs = if j == IDX_GLY {
                    &FLEX_DIHEDRAL_ANY_GLY
                } else if j == IDX_PRO {
                    if i == IDX_GLY {
                        &FLEX_DIHEDRAL_GLY_PRO
                    } else {
                        &FLEX_DIHEDRAL_OTHER_PRO
                    }
                } else {
                    &FLEX_DIHEDRAL_BY_AA1[i]
                };
                let _ = writeln!(s, "{prefix}{} {prefix}{} 1 7 \\", aa1.name, aa2.name);
                for &c in coeffs {
                    write_sci(&mut s, c, 15);
                }
                s.push_str("\n\n");
            }
        }
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_aa_handles_standard_and_aliases() {
        for aa in &AMINO_ACIDS {
            assert_eq!(canonical_aa(aa.name), Some(aa.name));
            assert_eq!(aa_params(aa.name).unwrap().name, aa.name);
        }
        for alias in ["HSD", "HSE", "HSP", "HID", "HIE", "HIP"] {
            assert_eq!(canonical_aa(alias), Some("HIS"));
        }
        assert_eq!(canonical_aa("CYM"), Some("CYS"));
        assert_eq!(canonical_aa("UNK"), Some("UNK"));
        assert_eq!(canonical_aa("POPC"), None);
        assert_eq!(canonical_aa("HOH"), None);
    }

    #[test]
    fn render_aicg2p_itp_contains_all_sections_and_reference_rows() {
        let itp = render_aicg2p_itp();
        assert!(itp.contains("[ defaults ]"));
        assert!(itp.contains("[ atomtypes ]"));
        assert!(itp.contains("[ flexible_local_angle ]"));
        assert!(itp.contains("[ flexible_local_dihedral_angle ]"));

        // Spot-check exact atomtypes, angle spline, and dihedral rows.
        assert!(itp.contains("  ALA  6  71.090   0.000     A     0.3943   2.5104"));
        assert!(itp.contains(" iMET  7 131.190   0.000     A     0.4424   2.5104"));
        assert!(itp.contains("ALA 1 10 \\\n     1.3090E+00     2.0920E+01     0.0000E+00 \\"));
        assert!(itp.contains(
            "GLY PRO 1 7 \\\n     1.0417E+01    -1.8882E+00     3.5531E+00    -1.6736E-02    -4.6652E-01     1.0887E+00     2.6610E-01"
        ));
        assert!(itp.contains(
            "iVAL iPRO 1 7 \\\n     1.1338E+01     2.0292E+00     3.2217E+00     3.9664E-01    -3.1271E+00     7.8450E-01    -1.5477E+00"
        ));
        // 40 angle blocks (20 AA + 20 iAA) and 800 dihedral blocks (400 + 400).
        assert_eq!(itp.matches(" 1 10 \\").count(), 40);
        assert_eq!(itp.matches(" 1 7 \\").count(), 800);
    }

    #[test]
    fn write_sci_formats_edge_values_without_allocating() {
        let fmt = |v: f64, w: usize| {
            let mut s = String::new();
            write_sci(&mut s, v, w);
            s
        };
        assert_eq!(fmt(0.0, 15), "     0.0000E+00");
        assert_eq!(fmt(-0.0, 15), "     0.0000E+00");
        assert_eq!(fmt(0.38017, 15), "     3.8017E-01");
        assert_eq!(fmt(-126.56, 15), "    -1.2656E+02");
        assert_eq!(fmt(1.5e-100, 15), "    1.5000E-100");
        assert_eq!(fmt(9.99999, 15), "     1.0000E+01");
    }
}
