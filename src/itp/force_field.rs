use serde::Deserialize;

#[derive(Debug, Default, Deserialize)]
pub struct Atom {
    pub name: String,
    pub e: f64,
    pub s: f64,
    pub w: Option<f64>,
    pub q: f64,
    pub m: f64,
    pub p: bool,
}

#[derive(Debug, Default, Deserialize)]
pub struct LipidAtom {
    pub name: String,
    pub id: u64,
}

#[derive(Debug, Default, Deserialize)]
pub struct LipidBond {
    pub ids: [u64; 2],
    pub k: f64,
    pub r0: f64,
}

#[derive(Debug, Default, Deserialize)]
pub struct LipidAngle {
    pub ids: [u64; 3],
    pub k: f64,
    pub t0: f64,
}

#[derive(Debug, Default, Deserialize)]
pub struct Lipid {
    pub name: String,
    pub atoms: Vec<LipidAtom>,
    pub bonds: Vec<LipidBond>,
    pub angles: Vec<LipidAngle>,
}

#[derive(Debug, Default, Deserialize)]
pub struct ForceField {
    pub atoms: Vec<Atom>,
    pub lipids: Vec<Lipid>,
}
