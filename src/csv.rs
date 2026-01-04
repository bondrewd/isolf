use serde::Deserialize;
use serde::de::DeserializeOwned;
use std::path::Path;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
pub struct Angle {
    pub lipid: String,
    pub bead1: String,
    pub bead2: String,
    pub bead3: String,
    pub t0: f64,
    pub k: f64,
}
#[derive(Debug, Deserialize)]
pub struct Bond {
    pub lipid: String,
    pub bead1: String,
    pub bead2: String,
    pub r0: f64,
    pub k: f64,
}

#[derive(Debug, Deserialize)]
pub struct Charge {
    pub bead: String,
    pub q: f64,
}

#[derive(Debug, Deserialize)]
pub struct Hp {
    pub bead: String,
    pub e: f64,
    pub s: f64,
    pub w: f64,
}

#[derive(Debug, Deserialize)]
pub struct Lj {
    pub bead: String,
    pub e: f64,
    pub s: f64,
}

#[derive(Debug, Deserialize)]
pub struct Mass {
    pub bead: String,
    pub m: f64,
}

#[derive(Debug, Deserialize)]
pub struct Wca {
    pub bead: String,
    pub e: f64,
    pub s: f64,
}

#[derive(Debug, Deserialize)]
pub struct Apl {
    pub lipid: String,
    pub mean30: f64,
    pub std30: f64,
    pub mean40: f64,
    pub std40: f64,
    pub mean50: f64,
    pub std50: f64,
}

#[derive(thiserror::Error, Debug)]
pub enum IsolfError {
    #[error("CSV error at {path}: {source}")]
    Csv { path: PathBuf, source: csv::Error },
}

impl IsolfError {
    pub fn csv(path: PathBuf, source: csv::Error) -> Self {
        IsolfError::Csv { path, source }
    }
}

#[derive(Debug)]
pub struct IsolfFf {
    pub angle: Vec<Angle>,
    pub bond: Vec<Bond>,
    pub charge: Vec<Charge>,
    pub hp: Vec<Hp>,
    pub lj: Vec<Lj>,
    pub mass: Vec<Mass>,
    pub wca: Vec<Wca>,
}

#[derive(Debug)]
pub struct IsolfCsv {
    pub ff: IsolfFf,
    pub apl: Vec<Apl>,
}

pub fn parse_csv<T, P>(path: P) -> Result<Vec<T>, IsolfError>
where
    T: DeserializeOwned,
    P: AsRef<Path>,
{
    let path = path.as_ref();
    csv::Reader::from_path(path)
        .map_err(|e| IsolfError::csv(path.into(), e))?
        .deserialize()
        .map(|res| res.map_err(|e| IsolfError::csv(path.into(), e)))
        .collect()
}

pub fn parse_isolf_csv() -> Result<IsolfCsv, IsolfError> {
    let root_path = env!("CARGO_MANIFEST_DIR");
    let data_path = PathBuf::from(root_path).join("data");

    Ok(IsolfCsv {
        ff: IsolfFf {
            angle: parse_csv(data_path.join("ff").join("angle.csv"))?,
            bond: parse_csv(data_path.join("ff").join("bond.csv"))?,
            charge: parse_csv(data_path.join("ff").join("charge.csv"))?,
            hp: parse_csv(data_path.join("ff").join("hp.csv"))?,
            lj: parse_csv(data_path.join("ff").join("lj.csv"))?,
            mass: parse_csv(data_path.join("ff").join("mass.csv"))?,
            wca: parse_csv(data_path.join("ff").join("wca.csv"))?,
        },
        apl: parse_csv(data_path.join("apl.csv"))?,
    })
}
