#[derive(Debug, Default, Clone)]
pub struct Output {
    pub path: String,
    pub period: u64,
}

impl Output {
    pub fn new(path: &str, period: u64) -> Self {
        Output {
            path: path.to_string(),
            period,
        }
    }
}
