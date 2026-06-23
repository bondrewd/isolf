//! Print the iSoLF `.itp` topology to stdout.
//!
//! Used to regenerate the golden fixture when a deliberate change alters the
//! topology: `cargo run --example generate_itp > tests/fixtures/isolf.itp`.

fn main() {
    print!("{}", isolf::force_field::ForceField::isolf().to_itp());
}
