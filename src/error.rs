//! Error types for building a membrane.

use std::fmt;

/// Something went wrong while parsing a composition or building a membrane.
#[derive(Debug, Clone, PartialEq)]
pub enum BuildError {
    /// A leaflet composition listed no lipids.
    EmptyComposition,

    /// A composition component had no lipid name (an empty entry, or a `=WEIGHT`
    /// with nothing before the `=`).
    MalformedComponent(String),

    /// A composition proportion could not be parsed (expected `a/b`, `x%`, or a
    /// number).
    InvalidProportion(String),

    /// Composition weights did not sum to a positive value.
    NonPositiveComposition,

    /// A lipid name was not found in the force field.
    UnknownLipid(String),

    /// A lipid referenced a bead type not present in the force field.
    UnknownBead(String),

    /// A lipid chain has no phosphate bead to orient it around.
    MissingPhosphate(String),

    /// A leaflet ended up with no lipids (e.g. a box too small for the chosen
    /// lipids).
    EmptyLeaflet,

    /// More lipids were requested than fit in the given box.
    BoxTooSmall { requested: usize, capacity: usize },

    /// The outer radius is too small to contain the bilayer's thickness.
    VesicleTooSmall { outer: f64, thickness: f64 },

    /// Both radii were given, but the gap between them can't hold the bilayer.
    VesicleTooThin {
        inner: f64,
        outer: f64,
        minimum: f64,
    },
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyComposition => f.write_str("composition is empty"),
            Self::MalformedComponent(s) => {
                write!(f, "composition component '{s}' has no lipid name")
            }
            Self::InvalidProportion(s) => write!(f, "invalid composition proportion '{s}'"),
            Self::NonPositiveComposition => {
                f.write_str("composition weights must sum to a positive value")
            }
            Self::UnknownLipid(s) => write!(f, "unknown lipid '{s}'"),
            Self::UnknownBead(s) => write!(f, "unknown bead '{s}'"),
            Self::MissingPhosphate(s) => write!(f, "lipid '{s}' has no phosphate bead"),
            Self::EmptyLeaflet => {
                f.write_str("leaflet would contain no lipids; increase the box size or lipid count")
            }
            Self::BoxTooSmall {
                requested,
                capacity,
            } => write!(
                f,
                "box holds at most {capacity} lipids per leaflet but {requested} were requested"
            ),
            Self::VesicleTooSmall { outer, thickness } => write!(
                f,
                "outer radius {outer:.2} nm is smaller than the bilayer thickness {thickness:.2} nm"
            ),
            Self::VesicleTooThin {
                inner,
                outer,
                minimum,
            } => write!(
                f,
                "the gap between the inner ({inner:.2} nm) and outer ({outer:.2} nm) radii is too thin \
                 for these lipids; the bilayer needs at least {minimum:.2} nm"
            ),
        }
    }
}

impl std::error::Error for BuildError {}
