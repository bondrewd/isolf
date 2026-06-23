//! Canonical iSoLF force-field parameter set.
//!
//! Generated from the reference `isolf.toml` parameter database; the numeric
//! values are the iSoLF v2.3.0 parameters. Edit the force field by
//! editing this module — it is the single source of truth.

// These are empirical force-field parameters; some happen to fall near math
// constants (e.g. a bond length of 0.5235 nm next to PI/6).
#![allow(clippy::approx_constant)]

use super::{BeadType, ForceField, HarmonicAngle, HarmonicBond, Interaction, Lipid};

impl ForceField {
    /// The canonical iSoLF (Implicit Solvent Lipid Force Field) parameter set,
    /// covering the standard phospholipid head groups and acyl chains.
    pub fn isolf() -> Self {
        Self {
            version: "2.3.0".to_string(),
            bead_types: bead_types(),
            lipids: lipids(),
        }
    }
}

fn bead_types() -> Vec<BeadType> {
    vec![
        BeadType {
            name: "PHA".into(),
            mass: 110.0066,
            charge: -1.0,
            interaction: Interaction::WeeksChandlerAndersen {
                epsilon: 1.913,
                sigma: 0.578,
            },
        },
        BeadType {
            name: "CHO".into(),
            mass: 87.166,
            charge: 1.0,
            interaction: Interaction::WeeksChandlerAndersen {
                epsilon: 1.75,
                sigma: 0.65,
            },
        },
        BeadType {
            name: "ETH".into(),
            mass: 45.0845,
            charge: 1.0,
            interaction: Interaction::WeeksChandlerAndersen {
                epsilon: 1.624,
                sigma: 0.438,
            },
        },
        BeadType {
            name: "GLC".into(),
            mass: 75.0878,
            charge: 0.0,
            interaction: Interaction::WeeksChandlerAndersen {
                epsilon: 1.612,
                sigma: 0.562,
            },
        },
        BeadType {
            name: "SRI".into(),
            mass: 88.0868,
            charge: 0.0,
            interaction: Interaction::LennardJones {
                epsilon: 1.656,
                sigma: 0.588,
            },
        },
        BeadType {
            name: "PHO".into(),
            mass: 108.9986,
            charge: -1.0,
            interaction: Interaction::WeeksChandlerAndersen {
                epsilon: 1.913,
                sigma: 0.578,
            },
        },
        BeadType {
            name: "MID".into(),
            mass: 143.1196,
            charge: 0.0,
            interaction: Interaction::WeeksChandlerAndersen {
                epsilon: 1.875,
                sigma: 0.507,
            },
        },
        BeadType {
            name: "DL1".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.233,
                sigma: 0.704,
                omega: 1.286,
            },
        },
        BeadType {
            name: "DL2".into(),
            mass: 142.286,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.233,
                sigma: 0.704,
                omega: 1.286,
            },
        },
        BeadType {
            name: "DM1".into(),
            mass: 168.324,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.307,
                sigma: 0.727,
                omega: 1.475,
            },
        },
        BeadType {
            name: "DM2".into(),
            mass: 170.34,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.307,
                sigma: 0.727,
                omega: 1.475,
            },
        },
        BeadType {
            name: "DP1".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.894,
                sigma: 0.71,
                omega: 1.406,
            },
        },
        BeadType {
            name: "DP2".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.894,
                sigma: 0.71,
                omega: 1.406,
            },
        },
        BeadType {
            name: "DP3".into(),
            mass: 114.232,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.894,
                sigma: 0.71,
                omega: 1.406,
            },
        },
        BeadType {
            name: "DO1".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.63,
                sigma: 0.726,
                omega: 1.265,
            },
        },
        BeadType {
            name: "DO2".into(),
            mass: 136.238,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.63,
                sigma: 0.726,
                omega: 1.265,
            },
        },
        BeadType {
            name: "DO3".into(),
            mass: 170.34,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.63,
                sigma: 0.726,
                omega: 1.265,
            },
        },
        BeadType {
            name: "DS1".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.366,
                sigma: 0.701,
                omega: 1.389,
            },
        },
        BeadType {
            name: "DS2".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.366,
                sigma: 0.701,
                omega: 1.389,
            },
        },
        BeadType {
            name: "DS3".into(),
            mass: 170.34,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 2.366,
                sigma: 0.701,
                omega: 1.389,
            },
        },
        BeadType {
            name: "PO1".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.719,
                sigma: 0.738,
                omega: 1.186,
            },
        },
        BeadType {
            name: "PO2".into(),
            mass: 138.254,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.719,
                sigma: 0.738,
                omega: 1.186,
            },
        },
        BeadType {
            name: "PO3".into(),
            mass: 142.286,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.719,
                sigma: 0.738,
                omega: 1.186,
            },
        },
        BeadType {
            name: "SO1".into(),
            mass: 140.27,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.806,
                sigma: 0.722,
                omega: 1.254,
            },
        },
        BeadType {
            name: "SO2".into(),
            mass: 138.254,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.806,
                sigma: 0.722,
                omega: 1.254,
            },
        },
        BeadType {
            name: "SO3".into(),
            mass: 170.34,
            charge: 0.0,
            interaction: Interaction::Isolf {
                epsilon: 1.806,
                sigma: 0.722,
                omega: 1.254,
            },
        },
    ]
}

fn lipids() -> Vec<Lipid> {
    vec![
        Lipid {
            name: "DLPA".into(),
            beads: vec!["PHA".into(), "MID".into(), "DL1".into(), "DL2".into()],
            bonds: vec![
                HarmonicBond {
                    length: 0.4685,
                    force_constant: 1763.4748,
                },
                HarmonicBond {
                    length: 0.4745,
                    force_constant: 2922.184,
                },
                HarmonicBond {
                    length: 0.5645,
                    force_constant: 2069.9689,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 173.6364,
                    force_constant: 7.1514,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 11.4009,
                },
            ],
        },
        Lipid {
            name: "DLPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "DL1".into(),
                "DL2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4025,
                    force_constant: 1984.7592,
                },
                HarmonicBond {
                    length: 0.4715,
                    force_constant: 1952.6966,
                },
                HarmonicBond {
                    length: 0.4675,
                    force_constant: 2471.3237,
                },
                HarmonicBond {
                    length: 0.5505,
                    force_constant: 1613.7962,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 140.9091,
                    force_constant: 3.7505,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.7101,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.4786,
                },
            ],
        },
        Lipid {
            name: "DLPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "DL1".into(),
                "DL2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3565,
                    force_constant: 4645.8799,
                },
                HarmonicBond {
                    length: 0.4735,
                    force_constant: 2031.9526,
                },
                HarmonicBond {
                    length: 0.4755,
                    force_constant: 2898.3888,
                },
                HarmonicBond {
                    length: 0.5665,
                    force_constant: 2154.0865,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 106.3636,
                    force_constant: 4.6831,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 7.4433,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 11.9031,
                },
            ],
        },
        Lipid {
            name: "DLPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "DL1".into(),
                "DL2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3815,
                    force_constant: 2556.8209,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 1909.5121,
                },
                HarmonicBond {
                    length: 0.4665,
                    force_constant: 2343.0142,
                },
                HarmonicBond {
                    length: 0.5455,
                    force_constant: 1441.4569,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 100.9091,
                    force_constant: 2.9059,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.4369,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 8.9696,
                },
            ],
        },
        Lipid {
            name: "DLPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "DL1".into(),
                "DL2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4425,
                    force_constant: 8885.3633,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 1831.8375,
                },
                HarmonicBond {
                    length: 0.4735,
                    force_constant: 2797.4588,
                },
                HarmonicBond {
                    length: 0.5635,
                    force_constant: 2094.5199,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 8.1513,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 7.3486,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 11.5875,
                },
            ],
        },
        Lipid {
            name: "DMPA".into(),
            beads: vec!["PHA".into(), "MID".into(), "DM1".into(), "DM2".into()],
            bonds: vec![
                HarmonicBond {
                    length: 0.4675,
                    force_constant: 1798.0765,
                },
                HarmonicBond {
                    length: 0.5265,
                    force_constant: 2211.4453,
                },
                HarmonicBond {
                    length: 0.6715,
                    force_constant: 1503.6403,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 6.7631,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.3727,
                },
            ],
        },
        Lipid {
            name: "DMPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "DM1".into(),
                "DM2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4035,
                    force_constant: 2029.6147,
                },
                HarmonicBond {
                    length: 0.4715,
                    force_constant: 1909.9275,
                },
                HarmonicBond {
                    length: 0.5235,
                    force_constant: 2050.6087,
                },
                HarmonicBond {
                    length: 0.6635,
                    force_constant: 1228.7037,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 140.9091,
                    force_constant: 3.7407,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.0334,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.1357,
                },
            ],
        },
        Lipid {
            name: "DMPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "DM1".into(),
                "DM2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3565,
                    force_constant: 4747.6185,
                },
                HarmonicBond {
                    length: 0.4745,
                    force_constant: 2197.485,
                },
                HarmonicBond {
                    length: 0.5295,
                    force_constant: 2411.0787,
                },
                HarmonicBond {
                    length: 0.6755,
                    force_constant: 1524.5372,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 106.3636,
                    force_constant: 4.8453,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.5227,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 13.0184,
                },
            ],
        },
        Lipid {
            name: "DMPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "DM1".into(),
                "DM2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3815,
                    force_constant: 2567.1502,
                },
                HarmonicBond {
                    length: 0.4715,
                    force_constant: 1987.6395,
                },
                HarmonicBond {
                    length: 0.5195,
                    force_constant: 1925.4607,
                },
                HarmonicBond {
                    length: 0.6515,
                    force_constant: 1029.6688,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 99.0909,
                    force_constant: 2.6501,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.6415,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.4798,
                },
            ],
        },
        Lipid {
            name: "DMPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "DM1".into(),
                "DM2".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4435,
                    force_constant: 9430.0853,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 1847.0322,
                },
                HarmonicBond {
                    length: 0.5305,
                    force_constant: 2297.7119,
                },
                HarmonicBond {
                    length: 0.6775,
                    force_constant: 1555.2015,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 8.1785,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.4593,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.7777,
                },
            ],
        },
        Lipid {
            name: "DPPA".into(),
            beads: vec![
                "PHA".into(),
                "MID".into(),
                "DP1".into(),
                "DP2".into(),
                "DP3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4725,
                    force_constant: 1810.4703,
                },
                HarmonicBond {
                    length: 0.4995,
                    force_constant: 4793.1979,
                },
                HarmonicBond {
                    length: 0.6345,
                    force_constant: 51782.5469,
                },
                HarmonicBond {
                    length: 0.5745,
                    force_constant: 49240.8982,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 171.8182,
                    force_constant: 9.4853,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 22.2413,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 21.1266,
                },
            ],
        },
        Lipid {
            name: "DPPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "DP1".into(),
                "DP2".into(),
                "DP3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4025,
                    force_constant: 2030.8916,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 1918.7328,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 2746.815,
                },
                HarmonicBond {
                    length: 0.5625,
                    force_constant: 2157.3367,
                },
                HarmonicBond {
                    length: 0.5005,
                    force_constant: 2101.9025,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 117.2727,
                    force_constant: 3.379,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.2353,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 11.5372,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.3884,
                },
            ],
        },
        Lipid {
            name: "DPPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "DP1".into(),
                "DP2".into(),
                "DP3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3555,
                    force_constant: 4597.7787,
                },
                HarmonicBond {
                    length: 0.4795,
                    force_constant: 2169.1758,
                },
                HarmonicBond {
                    length: 0.4985,
                    force_constant: 4478.8599,
                },
                HarmonicBond {
                    length: 0.6345,
                    force_constant: 51396.7936,
                },
                HarmonicBond {
                    length: 0.5745,
                    force_constant: 49529.0984,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 104.5455,
                    force_constant: 5.6466,
                },
                HarmonicAngle {
                    angle: 175.4545,
                    force_constant: 11.6392,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 30.5951,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 22.4752,
                },
            ],
        },
        Lipid {
            name: "DPPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "DP1".into(),
                "DP2".into(),
                "DP3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3805,
                    force_constant: 2575.8486,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 2097.6156,
                },
                HarmonicBond {
                    length: 0.4675,
                    force_constant: 2530.2906,
                },
                HarmonicBond {
                    length: 0.5555,
                    force_constant: 1822.9846,
                },
                HarmonicBond {
                    length: 0.4935,
                    force_constant: 1852.5981,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 102.7273,
                    force_constant: 3.1117,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 6.0854,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.5978,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 8.6826,
                },
            ],
        },
        Lipid {
            name: "DPPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "DP1".into(),
                "DP2".into(),
                "DP3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4425,
                    force_constant: 8992.8549,
                },
                HarmonicBond {
                    length: 0.4745,
                    force_constant: 2086.224,
                },
                HarmonicBond {
                    length: 0.4975,
                    force_constant: 4450.9854,
                },
                HarmonicBond {
                    length: 0.6335,
                    force_constant: 47005.4472,
                },
                HarmonicBond {
                    length: 0.5745,
                    force_constant: 49037.897,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 8.807,
                },
                HarmonicAngle {
                    angle: 171.8182,
                    force_constant: 10.1312,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 20.2075,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 19.9637,
                },
            ],
        },
        Lipid {
            name: "DOPA".into(),
            beads: vec![
                "PHA".into(),
                "MID".into(),
                "DO1".into(),
                "DO2".into(),
                "DO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4665,
                    force_constant: 1728.5872,
                },
                HarmonicBond {
                    length: 0.4645,
                    force_constant: 2309.8495,
                },
                HarmonicBond {
                    length: 0.5085,
                    force_constant: 1788.8395,
                },
                HarmonicBond {
                    length: 0.5575,
                    force_constant: 1157.2747,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.5459,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.6842,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.057,
                },
            ],
        },
        Lipid {
            name: "DOPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "DO1".into(),
                "DO2".into(),
                "DO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4025,
                    force_constant: 1990.9749,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 1892.745,
                },
                HarmonicBond {
                    length: 0.4615,
                    force_constant: 2095.9687,
                },
                HarmonicBond {
                    length: 0.5035,
                    force_constant: 1644.4215,
                },
                HarmonicBond {
                    length: 0.5475,
                    force_constant: 989.7797,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 140.9091,
                    force_constant: 3.7324,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.3352,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.759,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 4.5454,
                },
            ],
        },
        Lipid {
            name: "DOPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "DO1".into(),
                "DO2".into(),
                "DO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3565,
                    force_constant: 4675.3439,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 1890.4082,
                },
                HarmonicBond {
                    length: 0.4645,
                    force_constant: 2285.6905,
                },
                HarmonicBond {
                    length: 0.5155,
                    force_constant: 2036.3043,
                },
                HarmonicBond {
                    length: 0.5605,
                    force_constant: 1211.3682,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 108.1818,
                    force_constant: 4.5993,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.1163,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.9289,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 5.6215,
                },
            ],
        },
        Lipid {
            name: "DOPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "DO1".into(),
                "DO2".into(),
                "DO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3805,
                    force_constant: 2544.3856,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 1914.2659,
                },
                HarmonicBond {
                    length: 0.4575,
                    force_constant: 1927.2703,
                },
                HarmonicBond {
                    length: 0.4975,
                    force_constant: 1516.9694,
                },
                HarmonicBond {
                    length: 0.5375,
                    force_constant: 859.2079,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 100.9091,
                    force_constant: 3.0537,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 5.5036,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.275,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 4.2659,
                },
            ],
        },
        Lipid {
            name: "DOPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "DO1".into(),
                "DO2".into(),
                "DO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4425,
                    force_constant: 8829.7109,
                },
                HarmonicBond {
                    length: 0.4665,
                    force_constant: 1770.1901,
                },
                HarmonicBond {
                    length: 0.4655,
                    force_constant: 2335.7991,
                },
                HarmonicBond {
                    length: 0.5095,
                    force_constant: 1809.3645,
                },
                HarmonicBond {
                    length: 0.5575,
                    force_constant: 1165.3019,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 8.0045,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.3878,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.8988,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.1405,
                },
            ],
        },
        Lipid {
            name: "DSPA".into(),
            beads: vec![
                "PHA".into(),
                "MID".into(),
                "DS1".into(),
                "DS2".into(),
                "DS3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4725,
                    force_constant: 1841.0139,
                },
                HarmonicBond {
                    length: 0.4995,
                    force_constant: 4709.6117,
                },
                HarmonicBond {
                    length: 0.6345,
                    force_constant: 50774.0941,
                },
                HarmonicBond {
                    length: 0.7005,
                    force_constant: 39358.9782,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 175.4545,
                    force_constant: 7.7745,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 32.8992,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 70.0234,
                },
            ],
        },
        Lipid {
            name: "DSPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "DS1".into(),
                "DS2".into(),
                "DS3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4025,
                    force_constant: 1951.8821,
                },
                HarmonicBond {
                    length: 0.4765,
                    force_constant: 2029.1603,
                },
                HarmonicBond {
                    length: 0.4955,
                    force_constant: 4201.5328,
                },
                HarmonicBond {
                    length: 0.6335,
                    force_constant: 46812.9296,
                },
                HarmonicBond {
                    length: 0.7005,
                    force_constant: 39495.7646,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 140.9091,
                    force_constant: 4.1673,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.941,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 13.4998,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 15.8906,
                },
            ],
        },
        Lipid {
            name: "DSPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "DS1".into(),
                "DS2".into(),
                "DS3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3555,
                    force_constant: 4583.0591,
                },
                HarmonicBond {
                    length: 0.4805,
                    force_constant: 2205.1861,
                },
                HarmonicBond {
                    length: 0.4995,
                    force_constant: 4643.6881,
                },
                HarmonicBond {
                    length: 0.6345,
                    force_constant: 51429.2558,
                },
                HarmonicBond {
                    length: 0.7005,
                    force_constant: 40292.6788,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 106.3636,
                    force_constant: 6.0572,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.92,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 33.5578,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 61.3964,
                },
            ],
        },
        Lipid {
            name: "DSPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "DS1".into(),
                "DS2".into(),
                "DS3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3805,
                    force_constant: 2478.1795,
                },
                HarmonicBond {
                    length: 0.4735,
                    force_constant: 2004.6771,
                },
                HarmonicBond {
                    length: 0.4975,
                    force_constant: 4475.5104,
                },
                HarmonicBond {
                    length: 0.6345,
                    force_constant: 51166.9923,
                },
                HarmonicBond {
                    length: 0.7005,
                    force_constant: 39172.2941,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 133.6364,
                    force_constant: 6.4439,
                },
                HarmonicAngle {
                    angle: 175.4545,
                    force_constant: 7.3115,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 13.8743,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 16.5246,
                },
            ],
        },
        Lipid {
            name: "DSPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "DS1".into(),
                "DS2".into(),
                "DS3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4435,
                    force_constant: 9235.751,
                },
                HarmonicBond {
                    length: 0.4755,
                    force_constant: 2041.8932,
                },
                HarmonicBond {
                    length: 0.4975,
                    force_constant: 4255.9195,
                },
                HarmonicBond {
                    length: 0.6345,
                    force_constant: 51453.5621,
                },
                HarmonicBond {
                    length: 0.7005,
                    force_constant: 39994.3043,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 8.3988,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.5279,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 25.6652,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 48.8164,
                },
            ],
        },
        Lipid {
            name: "POPA".into(),
            beads: vec![
                "PHA".into(),
                "MID".into(),
                "PO1".into(),
                "PO2".into(),
                "PO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4665,
                    force_constant: 1711.8591,
                },
                HarmonicBond {
                    length: 0.4685,
                    force_constant: 2437.0221,
                },
                HarmonicBond {
                    length: 0.5415,
                    force_constant: 2056.8662,
                },
                HarmonicBond {
                    length: 0.5185,
                    force_constant: 481.4885,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 6.1382,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 11.9056,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.8275,
                },
            ],
        },
        Lipid {
            name: "POPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "PO1".into(),
                "PO2".into(),
                "PO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4025,
                    force_constant: 1994.8377,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 1871.2843,
                },
                HarmonicBond {
                    length: 0.4645,
                    force_constant: 2258.5233,
                },
                HarmonicBond {
                    length: 0.5355,
                    force_constant: 1928.5928,
                },
                HarmonicBond {
                    length: 0.5015,
                    force_constant: 383.4372,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 140.9091,
                    force_constant: 3.7422,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.5685,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.2469,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 4.8591,
                },
            ],
        },
        Lipid {
            name: "POPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "PO1".into(),
                "PO2".into(),
                "PO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3565,
                    force_constant: 4682.7105,
                },
                HarmonicBond {
                    length: 0.4715,
                    force_constant: 1945.7284,
                },
                HarmonicBond {
                    length: 0.4725,
                    force_constant: 2710.9747,
                },
                HarmonicBond {
                    length: 0.5475,
                    force_constant: 2345.5751,
                },
                HarmonicBond {
                    length: 0.5245,
                    force_constant: 524.6947,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 108.1818,
                    force_constant: 4.7933,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 6.9601,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.3173,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.3149,
                },
            ],
        },
        Lipid {
            name: "POPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "PO1".into(),
                "PO2".into(),
                "PO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3815,
                    force_constant: 2579.5595,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 1925.8955,
                },
                HarmonicBond {
                    length: 0.4635,
                    force_constant: 2226.4036,
                },
                HarmonicBond {
                    length: 0.5265,
                    force_constant: 1625.8103,
                },
                HarmonicBond {
                    length: 0.4945,
                    force_constant: 346.8038,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 99.0909,
                    force_constant: 2.8033,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.357,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.6056,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 4.421,
                },
            ],
        },
        Lipid {
            name: "POPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "PO1".into(),
                "PO2".into(),
                "PO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4425,
                    force_constant: 8900.3877,
                },
                HarmonicBond {
                    length: 0.4665,
                    force_constant: 1772.1281,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 2546.0923,
                },
                HarmonicBond {
                    length: 0.5415,
                    force_constant: 1994.4299,
                },
                HarmonicBond {
                    length: 0.5185,
                    force_constant: 495.7537,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 7.9485,
                },
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 7.0625,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.1286,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.8918,
                },
            ],
        },
        Lipid {
            name: "SOPA".into(),
            beads: vec![
                "PHA".into(),
                "MID".into(),
                "SO1".into(),
                "SO2".into(),
                "SO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4675,
                    force_constant: 1703.2501,
                },
                HarmonicBond {
                    length: 0.4695,
                    force_constant: 2431.6779,
                },
                HarmonicBond {
                    length: 0.5455,
                    force_constant: 2205.1838,
                },
                HarmonicBond {
                    length: 0.5885,
                    force_constant: 1486.6234,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 177.2727,
                    force_constant: 6.1522,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.2971,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.6138,
                },
            ],
        },
        Lipid {
            name: "SOPC".into(),
            beads: vec![
                "CHO".into(),
                "PHO".into(),
                "MID".into(),
                "SO1".into(),
                "SO2".into(),
                "SO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4025,
                    force_constant: 1982.4036,
                },
                HarmonicBond {
                    length: 0.4715,
                    force_constant: 1942.0011,
                },
                HarmonicBond {
                    length: 0.4645,
                    force_constant: 2257.9851,
                },
                HarmonicBond {
                    length: 0.5365,
                    force_constant: 1971.7266,
                },
                HarmonicBond {
                    length: 0.5745,
                    force_constant: 1223.3944,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 142.7273,
                    force_constant: 3.3095,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.566,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 10.4714,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.4456,
                },
            ],
        },
        Lipid {
            name: "SOPE".into(),
            beads: vec![
                "ETH".into(),
                "PHO".into(),
                "MID".into(),
                "SO1".into(),
                "SO2".into(),
                "SO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3565,
                    force_constant: 4651.5092,
                },
                HarmonicBond {
                    length: 0.4725,
                    force_constant: 1947.8874,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 2591.4312,
                },
                HarmonicBond {
                    length: 0.5485,
                    force_constant: 2413.7047,
                },
                HarmonicBond {
                    length: 0.5915,
                    force_constant: 1597.0852,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 106.3636,
                    force_constant: 4.5047,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.6445,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.792,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 8.1336,
                },
            ],
        },
        Lipid {
            name: "SOPG".into(),
            beads: vec![
                "GLC".into(),
                "PHO".into(),
                "MID".into(),
                "SO1".into(),
                "SO2".into(),
                "SO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.3815,
                    force_constant: 2610.6696,
                },
                HarmonicBond {
                    length: 0.4705,
                    force_constant: 1937.8007,
                },
                HarmonicBond {
                    length: 0.4635,
                    force_constant: 2214.6943,
                },
                HarmonicBond {
                    length: 0.5295,
                    force_constant: 1697.09,
                },
                HarmonicBond {
                    length: 0.5635,
                    force_constant: 1000.4587,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 100.9091,
                    force_constant: 2.9407,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 5.4076,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 9.8651,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.0359,
                },
            ],
        },
        Lipid {
            name: "SOPS".into(),
            beads: vec![
                "SRI".into(),
                "PHO".into(),
                "MID".into(),
                "SO1".into(),
                "SO2".into(),
                "SO3".into(),
            ],
            bonds: vec![
                HarmonicBond {
                    length: 0.4435,
                    force_constant: 9094.0328,
                },
                HarmonicBond {
                    length: 0.4675,
                    force_constant: 1810.0408,
                },
                HarmonicBond {
                    length: 0.4715,
                    force_constant: 2665.939,
                },
                HarmonicBond {
                    length: 0.5445,
                    force_constant: 2201.6502,
                },
                HarmonicBond {
                    length: 0.5845,
                    force_constant: 1361.3512,
                },
            ],
            angles: vec![
                HarmonicAngle {
                    angle: 128.1818,
                    force_constant: 7.7438,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 6.6411,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 12.4501,
                },
                HarmonicAngle {
                    angle: 179.0909,
                    force_constant: 7.599,
                },
            ],
        },
    ]
}
