#[derive(Debug, Default, Clone)]
pub struct Output {
    pub path: String,
    pub period: u64,
}

impl Output {
    #[must_use]
    pub fn new(path: &str, period: u64) -> Self {
        Self {
            path: path.to_string(),
            period,
        }
    }
}
