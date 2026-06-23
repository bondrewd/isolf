//! Strongly-typed representation of the iSoLF (Implicit Solvent Lipid Force
//! Field) coarse-grained model.
//!
//! The reference implementation stored the force field as a TOML file in which
//! each bead's properties were scattered across separate `mass`, `charge`,
//! `lj`, `wca` and `isolf` tables keyed by bead name. Building a topology then
//! required fragile string-keyed cross-lookups between those tables, and the
//! interaction class of a bead was only implied by *which* table it happened to
//! appear in.
//!
//! Here each bead instead carries a single [`Interaction`] enum value that holds
//! exactly the parameters its interaction class needs (so a tail bead's `omega`
//! cannot be set on a head bead, and a missing parameter is a compile error
//! rather than a runtime `KeyError`). Beads are stored once, in a single
//! canonical order, which removes the accidental section-dependent ordering of
//! the original TOML.

mod isolf;

/// Charges below this magnitude (in units of the elementary charge) are treated
/// as neutral when classifying bead pairs into interaction tables.
///
/// iSoLF bead charges are exactly `-1`, `0` or `+1`, so any threshold in
/// `(0, 1)` separates them; the explicit constant keeps the intent obvious.
pub const CHARGE_TOLERANCE: f64 = 0.01;

/// Ratio between a Lennard-Jones pair's cut-off distance and its `sigma`.
pub const LJ_CUTOFF_RATIO: f64 = 2.5;

/// 2^(1/6): an LJ/WCA pair potential reaches its minimum (and the WCA repulsion
/// vanishes) at `d = 2^(1/6)·σ_ij`. A bead's packing radius is therefore
/// `2^(1/6)·σ/2`, so two beads with arithmetic-mean σ mixing touch at the
/// potential minimum. This is the σ→radius rule used when packing lipids.
pub const POTENTIAL_MINIMUM_FACTOR: f64 = 1.122_462_048_309_372_8;

/// Non-bonded interaction model assigned to a coarse-grained bead.
///
/// The variant determines which parameter table of the topology a bead (and the
/// pairs it forms) is written to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Interaction {
    /// Attractive Lennard-Jones interaction, used by polar head beads (the
    /// phosphatidylserine head bead `SRI` in the current parameter set).
    LennardJones { epsilon: f64, sigma: f64 },
    /// Purely repulsive Weeks–Chandler–Andersen interaction, used by the
    /// remaining head and glycerol-linker beads.
    WeeksChandlerAndersen { epsilon: f64, sigma: f64 },
    /// iSoLF implicit-solvent attractive potential, used by acyl-tail beads.
    /// `omega` is the well-width parameter that encodes solvent-mediated
    /// cohesion and exists only for this variant.
    Isolf {
        epsilon: f64,
        sigma: f64,
        omega: f64,
    },
}

impl Interaction {
    /// Well depth `epsilon` (kJ/mol), common to every interaction class.
    pub fn epsilon(&self) -> f64 {
        match *self {
            Interaction::LennardJones { epsilon, .. }
            | Interaction::WeeksChandlerAndersen { epsilon, .. }
            | Interaction::Isolf { epsilon, .. } => epsilon,
        }
    }

    /// Particle diameter `sigma` (nm), common to every interaction class.
    pub fn sigma(&self) -> f64 {
        match *self {
            Interaction::LennardJones { sigma, .. }
            | Interaction::WeeksChandlerAndersen { sigma, .. }
            | Interaction::Isolf { sigma, .. } => sigma,
        }
    }
}

/// A coarse-grained bead type: one entry of the topology's `atomtypes` table.
#[derive(Debug, Clone, PartialEq)]
pub struct BeadType {
    /// Four-or-fewer-letter bead name, e.g. `PHO`, `MID`, `DL1`.
    pub name: String,
    /// Mass in g/mol.
    pub mass: f64,
    /// Charge in units of the elementary charge.
    pub charge: f64,
    /// Non-bonded interaction class and its parameters.
    pub interaction: Interaction,
}

impl BeadType {
    /// Whether the bead carries a non-negligible charge, used to decide whether
    /// a polar bead pairs with it attractively (LJ) or repulsively (WCA).
    pub fn is_charged(&self) -> bool {
        self.charge.abs() > CHARGE_TOLERANCE
    }
}

/// Harmonic bond between two consecutive beads of a lipid chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicBond {
    /// Equilibrium length `r0` (nm).
    pub length: f64,
    /// Force constant `k` (kJ·mol⁻¹·nm⁻²).
    pub force_constant: f64,
}

/// Harmonic angle between three consecutive beads of a lipid chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicAngle {
    /// Equilibrium angle `theta0` (degrees).
    pub angle: f64,
    /// Force constant `k` (kJ·mol⁻¹·rad⁻²).
    pub force_constant: f64,
}

/// Bead names that carry the lipid's phosphate group: `PHA` for phosphatidic
/// acid (which is itself the head bead) and `PHO` for the phosphate of every
/// other head group.
const PHOSPHATE_BEADS: [&str; 2] = ["PHA", "PHO"];

/// A coarse-grained lipid: a linear chain of beads with sequential connectivity.
///
/// The chain is bonded head-to-tail, so bond `i` joins beads `i` and `i + 1`
/// and angle `i` spans beads `i`, `i + 1` and `i + 2`. A well-formed lipid of
/// `n` beads therefore has `n - 1` bonds and `n - 2` angles; see
/// [`Lipid::is_well_formed`].
#[derive(Debug, Clone, PartialEq)]
pub struct Lipid {
    /// Four-letter lipid residue name, e.g. `POPC`.
    pub name: String,
    /// Bead-type names in chain order, from head group to tail terminus.
    pub beads: Vec<String>,
    /// Bonds along the chain, in chain order.
    pub bonds: Vec<HarmonicBond>,
    /// Angles along the chain, in chain order.
    pub angles: Vec<HarmonicAngle>,
}

impl Lipid {
    /// Whether the bond and angle counts are consistent with a linear chain of
    /// [`beads`](Lipid::beads).
    pub fn is_well_formed(&self) -> bool {
        let beads = self.beads.len();
        beads >= 2 && self.bonds.len() + 1 == beads && self.angles.len() + 2 == beads
    }

    /// Index (0-based) of the phosphate bead in the chain, around which the
    /// builder orients the lipid. `None` if the chain has no phosphate bead.
    pub fn phosphate_index(&self) -> Option<usize> {
        self.beads
            .iter()
            .position(|bead| PHOSPHATE_BEADS.contains(&bead.as_str()))
    }
}

/// A non-bonded pair entry with its Lorentz–Berthelot-combined parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedPair {
    /// First bead name.
    pub first: String,
    /// Second bead name.
    pub second: String,
    /// Combined well depth `epsilon` (kJ/mol).
    pub epsilon: f64,
    /// Combined diameter `sigma` (nm).
    pub sigma: f64,
}

/// Combine two beads' parameters with the Lorentz–Berthelot rules: the
/// geometric mean of the well depths and the arithmetic mean of the diameters.
pub fn lorentz_berthelot(a: &Interaction, b: &Interaction) -> (f64, f64) {
    let epsilon = (a.epsilon() * b.epsilon()).sqrt();
    let sigma = (a.sigma() + b.sigma()) * 0.5;
    (epsilon, sigma)
}

/// The complete iSoLF force field: bead types and the lipids built from them.
#[derive(Debug, Clone, PartialEq)]
pub struct ForceField {
    /// Force-field version string (used only for provenance in the output).
    pub version: String,
    /// Bead types in canonical order. This order is reused for every section of
    /// the generated topology.
    pub bead_types: Vec<BeadType>,
    /// Lipid definitions.
    pub lipids: Vec<Lipid>,
}

impl ForceField {
    /// Look up a bead type by name.
    pub fn bead_type(&self, name: &str) -> Option<&BeadType> {
        self.bead_types.iter().find(|bead| bead.name == name)
    }

    /// Look up a lipid by name.
    pub fn lipid(&self, name: &str) -> Option<&Lipid> {
        self.lipids.iter().find(|lipid| lipid.name == name)
    }

    /// Bead types interacting through the attractive Lennard-Jones potential.
    pub fn lennard_jones_beads(&self) -> impl Iterator<Item = &BeadType> {
        self.bead_types
            .iter()
            .filter(|bead| matches!(bead.interaction, Interaction::LennardJones { .. }))
    }

    /// Bead types interacting through the repulsive WCA potential.
    pub fn wca_beads(&self) -> impl Iterator<Item = &BeadType> {
        self.bead_types
            .iter()
            .filter(|bead| matches!(bead.interaction, Interaction::WeeksChandlerAndersen { .. }))
    }

    /// Bead types interacting through the iSoLF tail potential.
    pub fn isolf_beads(&self) -> impl Iterator<Item = &BeadType> {
        self.bead_types
            .iter()
            .filter(|bead| matches!(bead.interaction, Interaction::Isolf { .. }))
    }

    /// Attractive Lennard-Jones pairs of the topology's `cg_LJ_parameters`
    /// table: every polar–polar pair plus every polar–charged pair.
    ///
    /// Pairs are emitted in the canonical bead order, polar bead outermost, to
    /// match the layout produced by the reference generator.
    pub fn lennard_jones_pairs(&self) -> Vec<NonbondedPair> {
        let polar: Vec<&BeadType> = self.lennard_jones_beads().collect();
        let wca: Vec<&BeadType> = self.wca_beads().collect();

        let mut pairs = Vec::new();
        for a in &polar {
            for b in &polar {
                pairs.push(combine_pair(a, b));
            }
            for b in &wca {
                if b.is_charged() {
                    pairs.push(combine_pair(a, b));
                }
            }
        }
        pairs
    }

    /// Repulsive WCA pairs of the topology's `cg_WCA_parameters` table: every
    /// head–head and head–tail pair, plus the polar bead against the neutral
    /// head beads and against every tail bead.
    pub fn wca_pairs(&self) -> Vec<NonbondedPair> {
        let polar: Vec<&BeadType> = self.lennard_jones_beads().collect();
        let wca: Vec<&BeadType> = self.wca_beads().collect();
        let tails: Vec<&BeadType> = self.isolf_beads().collect();

        let mut pairs = Vec::new();
        for a in &wca {
            for b in &wca {
                pairs.push(combine_pair(a, b));
            }
            for b in &tails {
                pairs.push(combine_pair(a, b));
            }
        }
        for a in &polar {
            for b in &wca {
                if !b.is_charged() {
                    pairs.push(combine_pair(a, b));
                }
            }
            for b in &tails {
                pairs.push(combine_pair(a, b));
            }
        }
        pairs
    }
}

fn combine_pair(a: &BeadType, b: &BeadType) -> NonbondedPair {
    let (epsilon, sigma) = lorentz_berthelot(&a.interaction, &b.interaction);
    NonbondedPair {
        first: a.name.clone(),
        second: b.name.clone(),
        epsilon,
        sigma,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ff() -> ForceField {
        ForceField::isolf()
    }

    #[test]
    fn lorentz_berthelot_uses_geometric_and_arithmetic_means() {
        let a = Interaction::LennardJones {
            epsilon: 1.656,
            sigma: 0.588,
        };
        let b = Interaction::WeeksChandlerAndersen {
            epsilon: 1.913,
            sigma: 0.578,
        };
        let (epsilon, sigma) = lorentz_berthelot(&a, &b);
        assert!((epsilon - (1.656_f64 * 1.913).sqrt()).abs() < 1e-12);
        assert!((sigma - (0.588 + 0.578) * 0.5).abs() < 1e-12);
        // Matches the SRI–PHA row of the reference topology once rounded.
        assert_eq!(format!("{epsilon:>8.3} {sigma:>8.3}"), "   1.780    0.583");
    }

    #[test]
    fn bead_charge_classification() {
        let ff = ff();
        let charged: Vec<&str> = ff
            .bead_types
            .iter()
            .filter(|b| b.is_charged())
            .map(|b| b.name.as_str())
            .collect();
        assert_eq!(charged, ["PHA", "CHO", "ETH", "PHO"]);
    }

    #[test]
    fn interaction_class_counts() {
        let ff = ff();
        assert_eq!(ff.bead_types.len(), 26);
        assert_eq!(ff.lennard_jones_beads().count(), 1);
        assert_eq!(ff.wca_beads().count(), 6);
        assert_eq!(ff.isolf_beads().count(), 19);
    }

    #[test]
    fn pair_table_sizes() {
        let ff = ff();
        // polar×polar (1) + polar×charged-head (4)
        assert_eq!(ff.lennard_jones_pairs().len(), 5);
        // head×head (36) + head×tail (114) + polar×neutral-head (2) + polar×tail (19)
        assert_eq!(ff.wca_pairs().len(), 171);
    }

    #[test]
    fn every_lipid_is_a_well_formed_chain_of_known_beads() {
        let ff = ff();
        assert_eq!(ff.lipids.len(), 35);
        for lipid in &ff.lipids {
            assert!(lipid.is_well_formed(), "{} is not well formed", lipid.name);
            for bead in &lipid.beads {
                assert!(
                    ff.bead_type(bead).is_some(),
                    "{} references unknown bead {bead}",
                    lipid.name
                );
            }
        }
    }

    #[test]
    fn pairs_combine_underlying_parameters() {
        let ff = ff();
        let pha_dl1 = ff
            .wca_pairs()
            .into_iter()
            .find(|p| p.first == "PHA" && p.second == "DL1")
            .expect("PHA–DL1 pair");
        // PHA: WCA{1.913, 0.578}, DL1: Isolf{2.233, 0.704}
        assert_eq!(
            format!("{:>8.3} {:>8.3}", pha_dl1.epsilon, pha_dl1.sigma),
            "   2.067    0.641"
        );
    }
}
