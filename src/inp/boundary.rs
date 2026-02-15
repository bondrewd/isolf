#[derive(Debug, Default, Clone)]
pub struct BoxSize {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl BoxSize {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        BoxSize { x, y, z }
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
    pub fn nobc() -> Self {
        Boundary::NoBc
    }

    pub fn pbc() -> Self {
        Boundary::Pbc { box_size: None }
    }

    pub fn pbc_with_box_size(x: f64, y: f64, z: f64) -> Self {
        Boundary::Pbc {
            box_size: Some(BoxSize::new(x, y, z)),
        }
    }
}
