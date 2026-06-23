//! Leaflet compositions: the relative amounts of each lipid in a leaflet.
//!
//! A composition is parsed from per-lipid weights (`"1"`, `"1/2"`, `"50%"`, or a
//! plain number), normalized to proportions that sum to one, and then used to
//! apportion an integer number of lipids across the species.

use crate::error::BuildError;

/// One lipid species and its normalized proportion within a leaflet.
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    /// Lipid residue name (e.g. `POPC`), as it appears in the force field.
    pub name: String,
    /// Fraction of the leaflet made up of this lipid, in `[0, 1]`.
    pub proportion: f64,
}

/// A leaflet composition: lipid species with proportions summing to one,
/// kept in the order they were supplied.
#[derive(Debug, Clone, PartialEq)]
pub struct Composition {
    components: Vec<Component>,
}

/// Parse a single composition weight.
///
/// Accepts a ratio `"a/b"`, a percentage `"x%"`, or a plain decimal number. The
/// returned value is a relative weight; it is normalized later against the other
/// components of the leaflet.
///
/// # Errors
///
/// Returns [`BuildError::InvalidProportion`] if `spec` is not a non-negative,
/// finite ratio, percentage, or number.
pub fn parse_proportion(spec: &str) -> Result<f64, BuildError> {
    let spec = spec.trim();
    let invalid = || BuildError::InvalidProportion(spec.to_string());

    let value = if let Some((numerator, denominator)) = spec.split_once('/') {
        let numerator: f64 = numerator.trim().parse().map_err(|_| invalid())?;
        let denominator: f64 = denominator.trim().parse().map_err(|_| invalid())?;
        if denominator == 0.0 {
            return Err(invalid());
        }
        numerator / denominator
    } else if let Some(percent) = spec.strip_suffix('%') {
        percent.trim().parse::<f64>().map_err(|_| invalid())? / 100.0
    } else {
        spec.parse::<f64>().map_err(|_| invalid())?
    };

    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(invalid())
    }
}

impl Composition {
    /// Build a composition from `(lipid name, weight)` pairs, normalizing the
    /// weights so the proportions sum to one. Order is preserved.
    ///
    /// # Errors
    ///
    /// Returns [`BuildError::EmptyComposition`] if no weights are given, or
    /// [`BuildError::NonPositiveComposition`] if they do not sum to a positive
    /// value.
    pub fn from_weights(
        weights: impl IntoIterator<Item = (String, f64)>,
    ) -> Result<Self, BuildError> {
        let weights: Vec<(String, f64)> = weights.into_iter().collect();
        if weights.is_empty() {
            return Err(BuildError::EmptyComposition);
        }
        let total: f64 = weights.iter().map(|(_, weight)| weight).sum();
        if total <= 0.0 {
            return Err(BuildError::NonPositiveComposition);
        }
        let components = weights
            .into_iter()
            .map(|(name, weight)| Component {
                name,
                proportion: weight / total,
            })
            .collect();
        Ok(Self { components })
    }

    /// Parse a composition from a `NAME=WEIGHT,NAME=WEIGHT` specification, as
    /// supplied on the command line. A bare `NAME` (no `=WEIGHT`) defaults to a
    /// weight of 1. Lipid names are upper-cased and weights are parsed with
    /// [`parse_proportion`]; the result is normalized.
    ///
    /// # Errors
    ///
    /// Returns a [`BuildError`] if a component has an empty lipid name or a weight
    /// that does not parse.
    ///
    /// # Examples
    ///
    /// ```
    /// use isolf::composition::Composition;
    /// let c = Composition::parse("POPC=3,DOPS=1").unwrap();
    /// assert_eq!(c.components().len(), 2);
    /// // A bare name is the same as `NAME=1`.
    /// assert_eq!(Composition::parse("DOPC"), Composition::parse("DOPC=1"));
    /// ```
    pub fn parse(spec: &str) -> Result<Self, BuildError> {
        let mut weights = Vec::new();
        for component in spec.split(',') {
            // `NAME=WEIGHT` sets an explicit weight; a bare `NAME` defaults to 1.
            let (name, weight) = match component.split_once('=') {
                Some((name, weight)) => (name.trim(), parse_proportion(weight)?),
                None => (component.trim(), 1.0),
            };
            if name.is_empty() {
                return Err(BuildError::MalformedComponent(component.to_string()));
            }
            weights.push((name.to_uppercase(), weight));
        }
        Self::from_weights(weights)
    }

    /// The lipid species and their proportions.
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// Apportion `total` lipids across the species, returning a count per
    /// component (aligned with [`components`](Composition::components)) that sums
    /// to `total`.
    ///
    /// Each lipid in turn goes to the species whose running fraction sits
    /// furthest below its target proportion, so the counts track the requested
    /// ratios as closely as possible at every step. Ties favour the earlier
    /// component.
    pub fn partition(&self, total: usize) -> Vec<usize> {
        let species = self.components.len();
        let mut counts = vec![0usize; species];

        for placed in 1..=total {
            // Hand this lipid to the species whose fraction of the lipids placed
            // so far is furthest below its target proportion.
            let denominator = (placed - 1).max(1) as f64;
            let mut chosen = 0;
            let mut lowest = f64::INFINITY;
            for (k, component) in self.components.iter().enumerate() {
                let shortfall = counts[k] as f64 / denominator - component.proportion;
                if shortfall < lowest {
                    lowest = shortfall;
                    chosen = k;
                }
            }
            counts[chosen] += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ratios_percentages_and_numbers() {
        assert_eq!(parse_proportion("1").unwrap(), 1.0);
        assert_eq!(parse_proportion("1/2").unwrap(), 0.5);
        assert_eq!(parse_proportion("3/4").unwrap(), 0.75);
        assert_eq!(parse_proportion("50%").unwrap(), 0.5);
        assert_eq!(parse_proportion(" 0.25 ").unwrap(), 0.25);
    }

    #[test]
    fn rejects_malformed_or_negative_proportions() {
        assert!(parse_proportion("abc").is_err());
        assert!(parse_proportion("1/0").is_err());
        assert!(parse_proportion("-1").is_err());
        assert!(parse_proportion("1/").is_err());
    }

    #[test]
    fn parse_reads_name_weight_specs() {
        let composition = Composition::parse("popc=3,dops=1").unwrap();
        let components = composition.components();
        assert_eq!(components.len(), 2);
        // Names are upper-cased and weights normalized.
        assert_eq!(components[0].name, "POPC");
        assert_eq!(components[0].proportion, 0.75);
        assert_eq!(components[1].name, "DOPS");
        assert_eq!(components[1].proportion, 0.25);

        // A single component with a fractional weight is fine.
        let single = Composition::parse("DPPC=1/2").unwrap();
        assert_eq!(
            single.components(),
            [Component {
                name: "DPPC".to_string(),
                proportion: 1.0
            }]
        );
    }

    #[test]
    fn parse_defaults_a_bare_name_to_weight_one() {
        // A bare `NAME` is the same as `NAME=1`.
        assert_eq!(Composition::parse("DOPC"), Composition::parse("DOPC=1"));
        // A bare name mixes with weighted components at weight 1.
        let mix = Composition::parse("POPC,DOPC=3").unwrap();
        assert_eq!(mix.components()[0].proportion, 0.25); // 1 / (1 + 3)
        assert_eq!(mix.components()[1].proportion, 0.75);
        // An empty name is still rejected.
        assert!(Composition::parse("=1").is_err());
        assert!(Composition::parse("").is_err());
    }

    #[test]
    fn normalizes_weights_to_proportions() {
        let composition =
            Composition::from_weights([("DOPC".into(), 3.0), ("DOPS".into(), 1.0)]).unwrap();
        let proportions: Vec<f64> = composition
            .components()
            .iter()
            .map(|c| c.proportion)
            .collect();
        assert_eq!(proportions, [0.75, 0.25]);
    }

    #[test]
    fn empty_composition_is_rejected() {
        assert_eq!(
            Composition::from_weights([]).unwrap_err(),
            BuildError::EmptyComposition
        );
    }

    #[test]
    fn partition_sums_to_total_and_matches_proportions() {
        let composition =
            Composition::from_weights([("DOPC".into(), 1.0), ("DOPS".into(), 1.0)]).unwrap();
        let counts = composition.partition(100);
        assert_eq!(counts.iter().sum::<usize>(), 100);
        assert_eq!(counts, [50, 50]);

        let three =
            Composition::from_weights([("A".into(), 2.0), ("B".into(), 1.0), ("C".into(), 1.0)])
                .unwrap();
        let counts = three.partition(8);
        assert_eq!(counts.iter().sum::<usize>(), 8);
        assert_eq!(counts, [4, 2, 2]);
    }

    #[test]
    fn partition_of_pure_composition_is_all_one_species() {
        let composition = Composition::from_weights([("POPC".into(), 1.0)]).unwrap();
        assert_eq!(composition.partition(1024), [1024]);
    }
}
