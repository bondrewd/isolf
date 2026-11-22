//! Canonical iSoLF force-field parameter set.
//!
//! Numeric values are the iSoLF v2.3.0 parameters stored as static data tables
//! (keeping code size tiny) and materialized by [`ForceField::isolf`]. Edit the
//! force field by editing this module — it is the single source of truth.

// These are empirical force-field parameters; some happen to fall near math
// constants (e.g. a bond length of 0.5235 nm next to PI/6).
#![allow(clippy::approx_constant)]

use super::{BeadType, ForceField, HarmonicAngle, HarmonicBond, Interaction, Lipid};

struct BeadSpec {
    name: &'static str,
    mass: f64,
    charge: f64,
    interaction: Interaction,
}

struct LipidSpec {
    name: &'static str,
    beads: &'static [&'static str],
    bonds: &'static [(f64, f64)],
    angles: &'static [(f64, f64)],
}

impl ForceField {
    /// The canonical iSoLF (Implicit Solvent Lipid Force Field) parameter set,
    /// covering the standard phospholipid head groups and acyl chains.
    pub fn isolf() -> Self {
        Self {
            version: "2.3.0".to_string(),
            bead_types: BEAD_TYPES
                .iter()
                .map(|b| BeadType {
                    name: b.name.to_string(),
                    mass: b.mass,
                    charge: b.charge,
                    interaction: b.interaction,
                })
                .collect(),
            lipids: LIPIDS
                .iter()
                .map(|l| Lipid {
                    name: l.name.to_string(),
                    beads: l.beads.iter().map(|s| (*s).to_string()).collect(),
                    bonds: l
                        .bonds
                        .iter()
                        .map(|&(length, force_constant)| HarmonicBond {
                            length,
                            force_constant,
                        })
                        .collect(),
                    angles: l
                        .angles
                        .iter()
                        .map(|&(angle, force_constant)| HarmonicAngle {
                            angle,
                            force_constant,
                        })
                        .collect(),
                })
                .collect(),
        }
    }
}

const BEAD_TYPES: &[BeadSpec] = &[
    BeadSpec {
        name: "PHA",
        mass: 110.0066,
        charge: -1.0,
        interaction: Interaction::WeeksChandlerAndersen {
            epsilon: 1.913,
            sigma: 0.578,
        },
    },
    BeadSpec {
        name: "CHO",
        mass: 87.166,
        charge: 1.0,
        interaction: Interaction::WeeksChandlerAndersen {
            epsilon: 1.75,
            sigma: 0.65,
        },
    },
    BeadSpec {
        name: "ETH",
        mass: 45.0845,
        charge: 1.0,
        interaction: Interaction::WeeksChandlerAndersen {
            epsilon: 1.624,
            sigma: 0.438,
        },
    },
    BeadSpec {
        name: "GLC",
        mass: 75.0878,
        charge: 0.0,
        interaction: Interaction::WeeksChandlerAndersen {
            epsilon: 1.612,
            sigma: 0.562,
        },
    },
    BeadSpec {
        name: "SRI",
        mass: 88.0868,
        charge: 0.0,
        interaction: Interaction::LennardJones {
            epsilon: 1.656,
            sigma: 0.588,
        },
    },
    BeadSpec {
        name: "PHO",
        mass: 108.9986,
        charge: -1.0,
        interaction: Interaction::WeeksChandlerAndersen {
            epsilon: 1.913,
            sigma: 0.578,
        },
    },
    BeadSpec {
        name: "MID",
        mass: 143.1196,
        charge: 0.0,
        interaction: Interaction::WeeksChandlerAndersen {
            epsilon: 1.875,
            sigma: 0.507,
        },
    },
    BeadSpec {
        name: "DL1",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.233,
            sigma: 0.704,
            omega: 1.286,
        },
    },
    BeadSpec {
        name: "DL2",
        mass: 142.286,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.233,
            sigma: 0.704,
            omega: 1.286,
        },
    },
    BeadSpec {
        name: "DM1",
        mass: 168.324,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.307,
            sigma: 0.727,
            omega: 1.475,
        },
    },
    BeadSpec {
        name: "DM2",
        mass: 170.34,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.307,
            sigma: 0.727,
            omega: 1.475,
        },
    },
    BeadSpec {
        name: "DP1",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.894,
            sigma: 0.71,
            omega: 1.406,
        },
    },
    BeadSpec {
        name: "DP2",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.894,
            sigma: 0.71,
            omega: 1.406,
        },
    },
    BeadSpec {
        name: "DP3",
        mass: 114.232,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.894,
            sigma: 0.71,
            omega: 1.406,
        },
    },
    BeadSpec {
        name: "DO1",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.63,
            sigma: 0.726,
            omega: 1.265,
        },
    },
    BeadSpec {
        name: "DO2",
        mass: 136.238,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.63,
            sigma: 0.726,
            omega: 1.265,
        },
    },
    BeadSpec {
        name: "DO3",
        mass: 170.34,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.63,
            sigma: 0.726,
            omega: 1.265,
        },
    },
    BeadSpec {
        name: "DS1",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.366,
            sigma: 0.701,
            omega: 1.389,
        },
    },
    BeadSpec {
        name: "DS2",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.366,
            sigma: 0.701,
            omega: 1.389,
        },
    },
    BeadSpec {
        name: "DS3",
        mass: 170.34,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 2.366,
            sigma: 0.701,
            omega: 1.389,
        },
    },
    BeadSpec {
        name: "PO1",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.719,
            sigma: 0.738,
            omega: 1.186,
        },
    },
    BeadSpec {
        name: "PO2",
        mass: 138.254,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.719,
            sigma: 0.738,
            omega: 1.186,
        },
    },
    BeadSpec {
        name: "PO3",
        mass: 142.286,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.719,
            sigma: 0.738,
            omega: 1.186,
        },
    },
    BeadSpec {
        name: "SO1",
        mass: 140.27,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.806,
            sigma: 0.722,
            omega: 1.254,
        },
    },
    BeadSpec {
        name: "SO2",
        mass: 138.254,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.806,
            sigma: 0.722,
            omega: 1.254,
        },
    },
    BeadSpec {
        name: "SO3",
        mass: 170.34,
        charge: 0.0,
        interaction: Interaction::Isolf {
            epsilon: 1.806,
            sigma: 0.722,
            omega: 1.254,
        },
    },
];

const LIPIDS: &[LipidSpec] = &[
    LipidSpec {
        name: "DLPA",
        beads: &["PHA", "MID", "DL1", "DL2"],
        bonds: &[(0.4685, 1763.4748), (0.4745, 2922.184), (0.5645, 2069.9689)],
        angles: &[(173.6364, 7.1514), (179.0909, 11.4009)],
    },
    LipidSpec {
        name: "DLPC",
        beads: &["CHO", "PHO", "MID", "DL1", "DL2"],
        bonds: &[
            (0.4025, 1984.7592),
            (0.4715, 1952.6966),
            (0.4675, 2471.3237),
            (0.5505, 1613.7962),
        ],
        angles: &[(140.9091, 3.7505), (179.0909, 6.7101), (179.0909, 9.4786)],
    },
    LipidSpec {
        name: "DLPE",
        beads: &["ETH", "PHO", "MID", "DL1", "DL2"],
        bonds: &[
            (0.3565, 4645.8799),
            (0.4735, 2031.9526),
            (0.4755, 2898.3888),
            (0.5665, 2154.0865),
        ],
        angles: &[(106.3636, 4.6831), (177.2727, 7.4433), (179.0909, 11.9031)],
    },
    LipidSpec {
        name: "DLPG",
        beads: &["GLC", "PHO", "MID", "DL1", "DL2"],
        bonds: &[
            (0.3815, 2556.8209),
            (0.4695, 1909.5121),
            (0.4665, 2343.0142),
            (0.5455, 1441.4569),
        ],
        angles: &[(100.9091, 2.9059), (179.0909, 5.4369), (179.0909, 8.9696)],
    },
    LipidSpec {
        name: "DLPS",
        beads: &["SRI", "PHO", "MID", "DL1", "DL2"],
        bonds: &[
            (0.4425, 8885.3633),
            (0.4695, 1831.8375),
            (0.4735, 2797.4588),
            (0.5635, 2094.5199),
        ],
        angles: &[(128.1818, 8.1513), (177.2727, 7.3486), (179.0909, 11.5875)],
    },
    LipidSpec {
        name: "DMPA",
        beads: &["PHA", "MID", "DM1", "DM2"],
        bonds: &[
            (0.4675, 1798.0765),
            (0.5265, 2211.4453),
            (0.6715, 1503.6403),
        ],
        angles: &[(177.2727, 6.7631), (179.0909, 12.3727)],
    },
    LipidSpec {
        name: "DMPC",
        beads: &["CHO", "PHO", "MID", "DM1", "DM2"],
        bonds: &[
            (0.4035, 2029.6147),
            (0.4715, 1909.9275),
            (0.5235, 2050.6087),
            (0.6635, 1228.7037),
        ],
        angles: &[(140.9091, 3.7407), (179.0909, 7.0334), (179.0909, 10.1357)],
    },
    LipidSpec {
        name: "DMPE",
        beads: &["ETH", "PHO", "MID", "DM1", "DM2"],
        bonds: &[
            (0.3565, 4747.6185),
            (0.4745, 2197.485),
            (0.5295, 2411.0787),
            (0.6755, 1524.5372),
        ],
        angles: &[(106.3636, 4.8453), (179.0909, 7.5227), (179.0909, 13.0184)],
    },
    LipidSpec {
        name: "DMPG",
        beads: &["GLC", "PHO", "MID", "DM1", "DM2"],
        bonds: &[
            (0.3815, 2567.1502),
            (0.4715, 1987.6395),
            (0.5195, 1925.4607),
            (0.6515, 1029.6688),
        ],
        angles: &[(99.0909, 2.6501), (179.0909, 5.6415), (179.0909, 9.4798)],
    },
    LipidSpec {
        name: "DMPS",
        beads: &["SRI", "PHO", "MID", "DM1", "DM2"],
        bonds: &[
            (0.4435, 9430.0853),
            (0.4695, 1847.0322),
            (0.5305, 2297.7119),
            (0.6775, 1555.2015),
        ],
        angles: &[(128.1818, 8.1785), (179.0909, 7.4593), (179.0909, 12.7777)],
    },
    LipidSpec {
        name: "DPPA",
        beads: &["PHA", "MID", "DP1", "DP2", "DP3"],
        bonds: &[
            (0.4725, 1810.4703),
            (0.4995, 4793.1979),
            (0.6345, 51782.5469),
            (0.5745, 49240.8982),
        ],
        angles: &[(171.8182, 9.4853), (179.0909, 22.2413), (179.0909, 21.1266)],
    },
    LipidSpec {
        name: "DPPC",
        beads: &["CHO", "PHO", "MID", "DP1", "DP2", "DP3"],
        bonds: &[
            (0.4025, 2030.8916),
            (0.4705, 1918.7328),
            (0.4705, 2746.815),
            (0.5625, 2157.3367),
            (0.5005, 2101.9025),
        ],
        angles: &[
            (117.2727, 3.379),
            (179.0909, 7.2353),
            (179.0909, 11.5372),
            (179.0909, 9.3884),
        ],
    },
    LipidSpec {
        name: "DPPE",
        beads: &["ETH", "PHO", "MID", "DP1", "DP2", "DP3"],
        bonds: &[
            (0.3555, 4597.7787),
            (0.4795, 2169.1758),
            (0.4985, 4478.8599),
            (0.6345, 51396.7936),
            (0.5745, 49529.0984),
        ],
        angles: &[
            (104.5455, 5.6466),
            (175.4545, 11.6392),
            (179.0909, 30.5951),
            (179.0909, 22.4752),
        ],
    },
    LipidSpec {
        name: "DPPG",
        beads: &["GLC", "PHO", "MID", "DP1", "DP2", "DP3"],
        bonds: &[
            (0.3805, 2575.8486),
            (0.4705, 2097.6156),
            (0.4675, 2530.2906),
            (0.5555, 1822.9846),
            (0.4935, 1852.5981),
        ],
        angles: &[
            (102.7273, 3.1117),
            (177.2727, 6.0854),
            (179.0909, 10.5978),
            (179.0909, 8.6826),
        ],
    },
    LipidSpec {
        name: "DPPS",
        beads: &["SRI", "PHO", "MID", "DP1", "DP2", "DP3"],
        bonds: &[
            (0.4425, 8992.8549),
            (0.4745, 2086.224),
            (0.4975, 4450.9854),
            (0.6335, 47005.4472),
            (0.5745, 49037.897),
        ],
        angles: &[
            (128.1818, 8.807),
            (171.8182, 10.1312),
            (179.0909, 20.2075),
            (179.0909, 19.9637),
        ],
    },
    LipidSpec {
        name: "DOPA",
        beads: &["PHA", "MID", "DO1", "DO2", "DO3"],
        bonds: &[
            (0.4665, 1728.5872),
            (0.4645, 2309.8495),
            (0.5085, 1788.8395),
            (0.5575, 1157.2747),
        ],
        angles: &[(179.0909, 5.5459), (179.0909, 10.6842), (179.0909, 5.057)],
    },
    LipidSpec {
        name: "DOPC",
        beads: &["CHO", "PHO", "MID", "DO1", "DO2", "DO3"],
        bonds: &[
            (0.4025, 1990.9749),
            (0.4705, 1892.745),
            (0.4615, 2095.9687),
            (0.5035, 1644.4215),
            (0.5475, 989.7797),
        ],
        angles: &[
            (140.9091, 3.7324),
            (179.0909, 6.3352),
            (179.0909, 9.759),
            (179.0909, 4.5454),
        ],
    },
    LipidSpec {
        name: "DOPE",
        beads: &["ETH", "PHO", "MID", "DO1", "DO2", "DO3"],
        bonds: &[
            (0.3565, 4675.3439),
            (0.4695, 1890.4082),
            (0.4645, 2285.6905),
            (0.5155, 2036.3043),
            (0.5605, 1211.3682),
        ],
        angles: &[
            (108.1818, 4.5993),
            (179.0909, 6.1163),
            (179.0909, 10.9289),
            (177.2727, 5.6215),
        ],
    },
    LipidSpec {
        name: "DOPG",
        beads: &["GLC", "PHO", "MID", "DO1", "DO2", "DO3"],
        bonds: &[
            (0.3805, 2544.3856),
            (0.4695, 1914.2659),
            (0.4575, 1927.2703),
            (0.4975, 1516.9694),
            (0.5375, 859.2079),
        ],
        angles: &[
            (100.9091, 3.0537),
            (177.2727, 5.5036),
            (179.0909, 9.275),
            (179.0909, 4.2659),
        ],
    },
    LipidSpec {
        name: "DOPS",
        beads: &["SRI", "PHO", "MID", "DO1", "DO2", "DO3"],
        bonds: &[
            (0.4425, 8829.7109),
            (0.4665, 1770.1901),
            (0.4655, 2335.7991),
            (0.5095, 1809.3645),
            (0.5575, 1165.3019),
        ],
        angles: &[
            (128.1818, 8.0045),
            (179.0909, 6.3878),
            (179.0909, 10.8988),
            (179.0909, 5.1405),
        ],
    },
    LipidSpec {
        name: "DSPA",
        beads: &["PHA", "MID", "DS1", "DS2", "DS3"],
        bonds: &[
            (0.4725, 1841.0139),
            (0.4995, 4709.6117),
            (0.6345, 50774.0941),
            (0.7005, 39358.9782),
        ],
        angles: &[(175.4545, 7.7745), (179.0909, 32.8992), (179.0909, 70.0234)],
    },
    LipidSpec {
        name: "DSPC",
        beads: &["CHO", "PHO", "MID", "DS1", "DS2", "DS3"],
        bonds: &[
            (0.4025, 1951.8821),
            (0.4765, 2029.1603),
            (0.4955, 4201.5328),
            (0.6335, 46812.9296),
            (0.7005, 39495.7646),
        ],
        angles: &[
            (140.9091, 4.1673),
            (179.0909, 7.941),
            (179.0909, 13.4998),
            (179.0909, 15.8906),
        ],
    },
    LipidSpec {
        name: "DSPE",
        beads: &["ETH", "PHO", "MID", "DS1", "DS2", "DS3"],
        bonds: &[
            (0.3555, 4583.0591),
            (0.4805, 2205.1861),
            (0.4995, 4643.6881),
            (0.6345, 51429.2558),
            (0.7005, 40292.6788),
        ],
        angles: &[
            (106.3636, 6.0572),
            (179.0909, 9.92),
            (179.0909, 33.5578),
            (179.0909, 61.3964),
        ],
    },
    LipidSpec {
        name: "DSPG",
        beads: &["GLC", "PHO", "MID", "DS1", "DS2", "DS3"],
        bonds: &[
            (0.3805, 2478.1795),
            (0.4735, 2004.6771),
            (0.4975, 4475.5104),
            (0.6345, 51166.9923),
            (0.7005, 39172.2941),
        ],
        angles: &[
            (133.6364, 6.4439),
            (175.4545, 7.3115),
            (179.0909, 13.8743),
            (179.0909, 16.5246),
        ],
    },
    LipidSpec {
        name: "DSPS",
        beads: &["SRI", "PHO", "MID", "DS1", "DS2", "DS3"],
        bonds: &[
            (0.4435, 9235.751),
            (0.4755, 2041.8932),
            (0.4975, 4255.9195),
            (0.6345, 51453.5621),
            (0.7005, 39994.3043),
        ],
        angles: &[
            (128.1818, 8.3988),
            (179.0909, 7.5279),
            (179.0909, 25.6652),
            (179.0909, 48.8164),
        ],
    },
    LipidSpec {
        name: "POPA",
        beads: &["PHA", "MID", "PO1", "PO2", "PO3"],
        bonds: &[
            (0.4665, 1711.8591),
            (0.4685, 2437.0221),
            (0.5415, 2056.8662),
            (0.5185, 481.4885),
        ],
        angles: &[(177.2727, 6.1382), (179.0909, 11.9056), (179.0909, 5.8275)],
    },
    LipidSpec {
        name: "POPC",
        beads: &["CHO", "PHO", "MID", "PO1", "PO2", "PO3"],
        bonds: &[
            (0.4025, 1994.8377),
            (0.4705, 1871.2843),
            (0.4645, 2258.5233),
            (0.5355, 1928.5928),
            (0.5015, 383.4372),
        ],
        angles: &[
            (140.9091, 3.7422),
            (179.0909, 6.5685),
            (179.0909, 10.2469),
            (179.0909, 4.8591),
        ],
    },
    LipidSpec {
        name: "POPE",
        beads: &["ETH", "PHO", "MID", "PO1", "PO2", "PO3"],
        bonds: &[
            (0.3565, 4682.7105),
            (0.4715, 1945.7284),
            (0.4725, 2710.9747),
            (0.5475, 2345.5751),
            (0.5245, 524.6947),
        ],
        angles: &[
            (108.1818, 4.7933),
            (177.2727, 6.9601),
            (179.0909, 12.3173),
            (179.0909, 6.3149),
        ],
    },
    LipidSpec {
        name: "POPG",
        beads: &["GLC", "PHO", "MID", "PO1", "PO2", "PO3"],
        bonds: &[
            (0.3815, 2579.5595),
            (0.4705, 1925.8955),
            (0.4635, 2226.4036),
            (0.5265, 1625.8103),
            (0.4945, 346.8038),
        ],
        angles: &[
            (99.0909, 2.8033),
            (179.0909, 5.357),
            (179.0909, 9.6056),
            (179.0909, 4.421),
        ],
    },
    LipidSpec {
        name: "POPS",
        beads: &["SRI", "PHO", "MID", "PO1", "PO2", "PO3"],
        bonds: &[
            (0.4425, 8900.3877),
            (0.4665, 1772.1281),
            (0.4695, 2546.0923),
            (0.5415, 1994.4299),
            (0.5185, 495.7537),
        ],
        angles: &[
            (128.1818, 7.9485),
            (177.2727, 7.0625),
            (179.0909, 12.1286),
            (179.0909, 5.8918),
        ],
    },
    LipidSpec {
        name: "SOPA",
        beads: &["PHA", "MID", "SO1", "SO2", "SO3"],
        bonds: &[
            (0.4675, 1703.2501),
            (0.4695, 2431.6779),
            (0.5455, 2205.1838),
            (0.5885, 1486.6234),
        ],
        angles: &[(177.2727, 6.1522), (179.0909, 12.2971), (179.0909, 7.6138)],
    },
    LipidSpec {
        name: "SOPC",
        beads: &["CHO", "PHO", "MID", "SO1", "SO2", "SO3"],
        bonds: &[
            (0.4025, 1982.4036),
            (0.4715, 1942.0011),
            (0.4645, 2257.9851),
            (0.5365, 1971.7266),
            (0.5745, 1223.3944),
        ],
        angles: &[
            (142.7273, 3.3095),
            (179.0909, 6.566),
            (179.0909, 10.4714),
            (179.0909, 6.4456),
        ],
    },
    LipidSpec {
        name: "SOPE",
        beads: &["ETH", "PHO", "MID", "SO1", "SO2", "SO3"],
        bonds: &[
            (0.3565, 4651.5092),
            (0.4725, 1947.8874),
            (0.4705, 2591.4312),
            (0.5485, 2413.7047),
            (0.5915, 1597.0852),
        ],
        angles: &[
            (106.3636, 4.5047),
            (179.0909, 6.6445),
            (179.0909, 12.792),
            (179.0909, 8.1336),
        ],
    },
    LipidSpec {
        name: "SOPG",
        beads: &["GLC", "PHO", "MID", "SO1", "SO2", "SO3"],
        bonds: &[
            (0.3815, 2610.6696),
            (0.4705, 1937.8007),
            (0.4635, 2214.6943),
            (0.5295, 1697.09),
            (0.5635, 1000.4587),
        ],
        angles: &[
            (100.9091, 2.9407),
            (179.0909, 5.4076),
            (179.0909, 9.8651),
            (179.0909, 6.0359),
        ],
    },
    LipidSpec {
        name: "SOPS",
        beads: &["SRI", "PHO", "MID", "SO1", "SO2", "SO3"],
        bonds: &[
            (0.4435, 9094.0328),
            (0.4675, 1810.0408),
            (0.4715, 2665.939),
            (0.5445, 2201.6502),
            (0.5845, 1361.3512),
        ],
        angles: &[
            (128.1818, 7.7438),
            (179.0909, 6.6411),
            (179.0909, 12.4501),
            (179.0909, 7.599),
        ],
    },
];
