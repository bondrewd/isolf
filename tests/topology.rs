//! End-to-end check that the rendered iSoLF topology stays byte-for-byte stable.
//!
//! The fixture in `tests/fixtures/isolf.itp` was produced by this crate and
//! verified to be value-identical to the reference `gen_itp.py` output (it
//! differs only in the version string and in using a single consistent bead
//! ordering). This test guards against accidental regressions in the renderer
//! or the parameter data.

use isolf::force_field::ForceField;

#[test]
fn matches_golden_snapshot() {
    let generated = ForceField::isolf().to_itp();
    let expected = include_str!("fixtures/isolf.itp");
    assert_eq!(
        generated, expected,
        "rendered topology drifted from tests/fixtures/isolf.itp; \
         regenerate with `cargo run --example generate_itp` if the change is intended"
    );
}
