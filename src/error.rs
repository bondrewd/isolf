//! Error types for building a membrane.

/// Something went wrong while parsing a composition or building a membrane.
#[derive(thiserror::Error, Debug, Clone, PartialEq)]
pub enum BuildError {
    /// A leaflet composition listed no lipids.
    #[error("composition is empty")]
    EmptyComposition,

    /// A composition component had no lipid name (an empty entry, or a `=WEIGHT`
    /// with nothing before the `=`).
    #[error("composition component '{0}' has no lipid name")]
    MalformedComponent(String),

    /// A composition proportion could not be parsed (expected `a/b`, `x%`, or a
    /// number).
    #[error("invalid composition proportion '{0}'")]
    InvalidProportion(String),

    /// Composition weights did not sum to a positive value.
    #[error("composition weights must sum to a positive value")]
    NonPositiveComposition,

    /// A lipid name was not found in the force field.
    #[error("unknown lipid '{0}'")]
    UnknownLipid(String),

    /// A lipid referenced a bead type not present in the force field.
    #[error("unknown bead '{0}'")]
    UnknownBead(String),

    /// A lipid chain has no phosphate bead to orient it around.
    #[error("lipid '{0}' has no phosphate bead")]
    MissingPhosphate(String),

    /// A leaflet ended up with no lipids (e.g. a box too small for the chosen
    /// lipids).
    #[error("leaflet would contain no lipids; increase the box size or lipid count")]
    EmptyLeaflet,

    /// More lipids were requested than fit in the given box.
    #[error("box holds at most {capacity} lipids per leaflet but {requested} were requested")]
    BoxTooSmall { requested: usize, capacity: usize },

    /// The outer radius is too small to contain the bilayer's thickness.
    #[error("outer radius {outer:.2} nm is smaller than the bilayer thickness {thickness:.2} nm")]
    VesicleTooSmall { outer: f64, thickness: f64 },

    /// Both radii were given, but the gap between them can't hold the bilayer.
    #[error(
        "the gap between the inner ({inner:.2} nm) and outer ({outer:.2} nm) radii is too thin \
         for these lipids; the bilayer needs at least {minimum:.2} nm"
    )]
    VesicleTooThin {
        inner: f64,
        outer: f64,
        minimum: f64,
    },
}
