//! GENESIS topology (`.top`) rendering of a [`Membrane`].
//!
//! Lists the molecule counts (in the shared canonical order, matching the
//! `.gro`), `#include`s the iSoLF force-field topology, and enables
//! electrostatics across every lipid chain.

use std::fmt;

use crate::membrane::Membrane;

/// Wraps a [`Membrane`] for `.top` rendering via [`fmt::Display`].
pub struct Top<'a> {
    membrane: &'a Membrane,
}

impl<'a> Top<'a> {
    /// Wrap a membrane for rendering.
    pub fn new(membrane: &'a Membrane) -> Self {
        Self { membrane }
    }
}

impl Membrane {
    /// Render this membrane as a GENESIS `.top` topology file.
    pub fn to_top(&self) -> String {
        Top::new(self).to_string()
    }
}

impl fmt::Display for Top<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let membrane = self.membrane;
        writeln!(f, "; Implicit Solvent Lipid Forcefield (iSoLF)")?;
        writeln!(f, "#include \"./isolf.itp\"")?;
        writeln!(f)?;

        writeln!(f, "[ system ]")?;
        writeln!(f, "{}", membrane.name)?;
        writeln!(f)?;

        writeln!(f, "[ molecules ]")?;
        for (lipid, count) in &membrane.lipid_counts {
            writeln!(f, "{lipid} {count}")?;
        }
        writeln!(f)?;

        writeln!(f, "[ cg_ele_chain_pairs ]")?;
        let total = membrane.total_lipids();
        writeln!(f, "ON 1 - {total} : 1 - {total}")
    }
}
