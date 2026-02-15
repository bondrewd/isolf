#[derive(thiserror::Error, Debug)]
pub enum IsolfError {
    #[error("Invalid input '{input}': {reason}")]
    InvalidInput { input: String, reason: String },

    #[error("Missing required field '{field}'")]
    MissingRequiredField { field: String },

    #[error("Invalid value '{value}' for field '{field}': {reason}")]
    InvalidFieldValue {
        field: String,
        value: String,
        reason: String,
    },
}

impl IsolfError {
    pub fn invalid_input(input: String, reason: String) -> Self {
        IsolfError::InvalidInput { input, reason }
    }

    pub fn missing_required_field<F: ToString>(field: F) -> Self {
        IsolfError::MissingRequiredField {
            field: field.to_string(),
        }
    }

    pub fn invalid_field_value<F: ToString, V: ToString, R: ToString>(
        field: F,
        value: V,
        reason: R,
    ) -> Self {
        IsolfError::InvalidFieldValue {
            field: field.to_string(),
            value: value.to_string(),
            reason: reason.to_string(),
        }
    }
}
