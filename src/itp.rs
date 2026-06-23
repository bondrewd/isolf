//! Rendering of a [`ForceField`] to a GROMACS-style iSoLF topology (`.itp`).
//!
//! The output is consumed by GENESIS and is column-sensitive, so the layout
//! here mirrors the fixed-width formatting of the reference generator. The file
//! is built from the force field in five parts:
//!
//! 1. `[ atomtypes ]` — one row per bead type.
//! 2. `[ cg_LJ_parameters ]` — attractive Lennard-Jones pairs (with cut-off).
//! 3. `[ cg_WCA_parameters ]` — repulsive WCA pairs.
//! 4. `[ cg_ISOLF_parameters ]` — per-bead iSoLF tail parameters.
//! 5. one `[ moleculetype ]` block per lipid (`atoms`, `bonds`, `angles`).

use std::fmt;

use crate::force_field::{ForceField, Interaction, LJ_CUTOFF_RATIO, Lipid};

/// Number of bonded neighbours excluded from non-bonded interactions
/// (`nrexcl`); iSoLF excludes up to 1-3 neighbours.
const EXCLUDED_NEIGHBORS: u32 = 2;

/// A [`ForceField`] paired with a [`fmt::Display`] implementation that renders
/// it as an iSoLF `.itp` topology.
///
/// ```
/// use isolf::force_field::ForceField;
/// use isolf::itp::Itp;
///
/// let topology = Itp::new(&ForceField::isolf()).to_string();
/// assert!(topology.starts_with("; "));
/// ```
pub struct Itp<'a> {
    force_field: &'a ForceField,
}

impl<'a> Itp<'a> {
    /// Wrap a force field for rendering.
    pub fn new(force_field: &'a ForceField) -> Self {
        Self { force_field }
    }
}

impl ForceField {
    /// Render this force field as an iSoLF `.itp` topology string.
    pub fn to_itp(&self) -> String {
        Itp::new(self).to_string()
    }
}

impl fmt::Display for Itp<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ff = self.force_field;
        write_header(f, &ff.version)?;
        write_atomtypes(f, ff)?;
        write_lennard_jones(f, ff)?;
        write_wca(f, ff)?;
        write_isolf(f, ff)?;
        for lipid in &ff.lipids {
            write_lipid(f, ff, lipid)?;
        }
        Ok(())
    }
}

fn write_header(f: &mut fmt::Formatter<'_>, version: &str) -> fmt::Result {
    writeln!(f, "; ------------------------------------------------")?;
    writeln!(
        f,
        "; ISOLF: Implicit Solvent Lipid Force Field v{version:>5}"
    )?;
    writeln!(f, "; ------------------------------------------------")?;
    writeln!(f)
}

fn write_atomtypes(f: &mut fmt::Formatter<'_>, ff: &ForceField) -> fmt::Result {
    writeln!(f, "[ atomtypes ]")?;
    writeln!(f, "; name   n     mass   charge ptype      rmin      eps")?;
    writeln!(f, ";    -   -    g/mol        e     -        nm   kJ/mol")?;
    for bead in &ff.bead_types {
        writeln!(
            f,
            "{:>6}   1 {:>8.4} {:>8.3}     A    0.0000   0.0000",
            bead.name, bead.mass, bead.charge
        )?;
    }
    writeln!(f)
}

fn write_lennard_jones(f: &mut fmt::Formatter<'_>, ff: &ForceField) -> fmt::Result {
    writeln!(f, "[ cg_LJ_parameters ]")?;
    writeln!(f, "; name  name  epsilon    sigma  cut-off")?;
    writeln!(f, ";    -     -   kJ/mol       nm       nm")?;
    for pair in ff.lennard_jones_pairs() {
        let cutoff = pair.sigma * LJ_CUTOFF_RATIO;
        writeln!(
            f,
            "{:>6}{:>6} {:>8.3} {:>8.3} {:>8.3}",
            pair.first, pair.second, pair.epsilon, pair.sigma, cutoff
        )?;
    }
    writeln!(f)
}

fn write_wca(f: &mut fmt::Formatter<'_>, ff: &ForceField) -> fmt::Result {
    writeln!(f, "[ cg_WCA_parameters ]")?;
    writeln!(f, "; name  name  epsilon    sigma")?;
    writeln!(f, ";    -     -   kJ/mol       nm")?;
    for pair in ff.wca_pairs() {
        writeln!(
            f,
            "{:>6}{:>6} {:>8.3} {:>8.3}",
            pair.first, pair.second, pair.epsilon, pair.sigma
        )?;
    }
    writeln!(f)
}

fn write_isolf(f: &mut fmt::Formatter<'_>, ff: &ForceField) -> fmt::Result {
    writeln!(f, "[ cg_ISOLF_parameters ]")?;
    writeln!(f, "; name  epsilon    sigma    omega")?;
    writeln!(f, ";    -   kJ/mol       nm       nm")?;
    for bead in ff.isolf_beads() {
        if let Interaction::Isolf {
            epsilon,
            sigma,
            omega,
        } = bead.interaction
        {
            writeln!(
                f,
                "{:>6} {:>8.3} {:>8.3} {:>8.3}",
                bead.name, epsilon, sigma, omega
            )?;
        }
    }
    writeln!(f)
}

fn write_lipid(f: &mut fmt::Formatter<'_>, ff: &ForceField, lipid: &Lipid) -> fmt::Result {
    writeln!(f, "[ moleculetype ]")?;
    writeln!(f, "; name  nrexcl")?;
    writeln!(f, "{:>6}       {}", lipid.name, EXCLUDED_NEIGHBORS)?;
    writeln!(f)?;

    writeln!(f, "[ atoms ]")?;
    writeln!(f, "; nr   type  resnr    res   atom   cg   charge     mass")?;
    writeln!(f, ";  -      -      -      -      -    -        e      amu")?;
    for (index, bead_name) in lipid.beads.iter().enumerate() {
        // Validated by `Lipid` construction / `ForceField::isolf`; an unknown
        // bead here means a malformed force field, so surface it loudly.
        let bead = ff
            .bead_type(bead_name)
            .unwrap_or_else(|| panic!("lipid {} references unknown bead {bead_name}", lipid.name));
        writeln!(
            f,
            "{:>4} {:>6}      1 {:>6} {:>6}    1 {:>8.3} {:>8.4}",
            index + 1,
            bead.name,
            lipid.name,
            bead.name,
            bead.charge,
            bead.mass
        )?;
    }
    writeln!(f)?;

    writeln!(f, "[ bonds ]")?;
    writeln!(f, ";  i   j   f       eq           coef")?;
    writeln!(f, ";  -   -   -       nm  kJ*nm-2*mol-1")?;
    for (index, bond) in lipid.bonds.iter().enumerate() {
        writeln!(
            f,
            "{:>4}{:>4}   1 {:>8.4} {:>14.4}",
            index + 1,
            index + 2,
            bond.length,
            bond.force_constant
        )?;
    }
    writeln!(f)?;

    writeln!(f, "[ angles ]")?;
    writeln!(f, ";  i   j   k   f       eq           coef")?;
    writeln!(f, ";  -   -   -   -      deg kJ*rad-2*mol-1")?;
    for (index, angle) in lipid.angles.iter().enumerate() {
        writeln!(
            f,
            "{:>4}{:>4}{:>4}   1 {:>8.4} {:>14.4}",
            index + 1,
            index + 2,
            index + 3,
            angle.angle,
            angle.force_constant
        )?;
    }
    writeln!(f)
}

#[cfg(test)]
mod tests {
    use crate::force_field::ForceField;

    fn topology() -> String {
        ForceField::isolf().to_itp()
    }

    #[test]
    fn renders_the_full_atomtypes_block() {
        // Consistent canonical bead order, so this block matches the reference
        // generator byte-for-byte.
        let expected = "\
[ atomtypes ]
; name   n     mass   charge ptype      rmin      eps
;    -   -    g/mol        e     -        nm   kJ/mol
   PHA   1 110.0066   -1.000     A    0.0000   0.0000
   CHO   1  87.1660    1.000     A    0.0000   0.0000
   ETH   1  45.0845    1.000     A    0.0000   0.0000
   GLC   1  75.0878    0.000     A    0.0000   0.0000
   SRI   1  88.0868    0.000     A    0.0000   0.0000
   PHO   1 108.9986   -1.000     A    0.0000   0.0000
   MID   1 143.1196    0.000     A    0.0000   0.0000
";
        assert!(topology().contains(expected), "atomtypes block mismatch");
    }

    #[test]
    fn renders_the_full_lennard_jones_block() {
        // The LJ table involves no tail beads, so its order matches the
        // reference exactly.
        let expected = "\
[ cg_LJ_parameters ]
; name  name  epsilon    sigma  cut-off
;    -     -   kJ/mol       nm       nm
   SRI   SRI    1.656    0.588    1.470
   SRI   PHA    1.780    0.583    1.458
   SRI   CHO    1.702    0.619    1.547
   SRI   ETH    1.640    0.513    1.282
   SRI   PHO    1.780    0.583    1.458
";
        assert!(topology().contains(expected), "LJ block mismatch");
    }

    #[test]
    fn renders_dlpc_moleculetype_block() {
        let expected = "\
[ moleculetype ]
; name  nrexcl
  DLPC       2

[ atoms ]
; nr   type  resnr    res   atom   cg   charge     mass
;  -      -      -      -      -    -        e      amu
   1    CHO      1   DLPC    CHO    1    1.000  87.1660
   2    PHO      1   DLPC    PHO    1   -1.000 108.9986
   3    MID      1   DLPC    MID    1    0.000 143.1196
   4    DL1      1   DLPC    DL1    1    0.000 140.2700
   5    DL2      1   DLPC    DL2    1    0.000 142.2860

[ bonds ]
;  i   j   f       eq           coef
;  -   -   -       nm  kJ*nm-2*mol-1
   1   2   1   0.4025      1984.7592
   2   3   1   0.4715      1952.6966
   3   4   1   0.4675      2471.3237
   4   5   1   0.5505      1613.7962

[ angles ]
;  i   j   k   f       eq           coef
;  -   -   -   -      deg kJ*rad-2*mol-1
   1   2   3   1 140.9091         3.7505
   2   3   4   1 179.0909         6.7101
   3   4   5   1 179.0909         9.4786
";
        assert!(topology().contains(expected), "DLPC block mismatch");
    }

    #[test]
    fn renders_dppc_block_with_three_tail_beads() {
        // Exercises a six-bead lipid (five bonds, four angles) and the
        // five-significant-figure force constants of the palmitoyl tail.
        let expected = "\
   1    CHO      1   DPPC    CHO    1    1.000  87.1660
   2    PHO      1   DPPC    PHO    1   -1.000 108.9986
   3    MID      1   DPPC    MID    1    0.000 143.1196
   4    DP1      1   DPPC    DP1    1    0.000 140.2700
   5    DP2      1   DPPC    DP2    1    0.000 140.2700
   6    DP3      1   DPPC    DP3    1    0.000 114.2320

[ bonds ]
;  i   j   f       eq           coef
;  -   -   -       nm  kJ*nm-2*mol-1
   1   2   1   0.4025      2030.8916
   2   3   1   0.4705      1918.7328
   3   4   1   0.4705      2746.8150
   4   5   1   0.5625      2157.3367
   5   6   1   0.5005      2101.9025
";
        assert!(topology().contains(expected), "DPPC block mismatch");
    }

    #[test]
    fn contains_reference_wca_and_isolf_values() {
        let topology = topology();
        // Spot-check representative rows whose values are fixed by the
        // reference, independent of ordering.
        for line in [
            "   PHA   DL1    2.067    0.641",    // head–tail WCA
            "   MID   MID    1.875    0.507",    // head–head WCA (self)
            "   SRI   GLC    1.634    0.575",    // polar–neutral-head WCA
            "   SRI   DL1    1.923    0.646",    // polar–tail WCA
            "   DP1    1.894    0.710    1.406", // iSoLF tail entry
            "   SO3    1.806    0.722    1.254",
        ] {
            assert!(topology.contains(line), "missing reference line: {line:?}");
        }
    }

    #[test]
    fn has_one_moleculetype_per_lipid() {
        assert_eq!(topology().matches("[ moleculetype ]").count(), 35);
    }
}
