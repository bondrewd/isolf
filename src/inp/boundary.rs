#[derive(Debug, Default, Clone)]
pub struct BoxSize {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl BoxSize {
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

#[derive(Debug, Default, Clone)]
pub enum Boundary {
    #[default]
    NoBc,
    Pbc {
        box_size: Option<BoxSize>,
    },
}

impl Boundary {
    #[must_use]
    pub const fn nobc() -> Self {
        Self::NoBc
    }

    #[must_use]
    pub const fn pbc() -> Self {
        Self::Pbc { box_size: None }
    }

    #[must_use]
    pub const fn pbc_with_box_size(x: f64, y: f64, z: f64) -> Self {
        Self::Pbc {
            box_size: Some(BoxSize::new(x, y, z)),
        }
    }
}
