//! `isolf lipids` table printer and Levenshtein lipid-name suggestions.

use std::error::Error;

use isolf::error::BuildError;
use isolf::force_field::ForceField;

/// Turn an [`UnknownLipid`](BuildError::UnknownLipid) into an actionable message
/// with close matches and a pointer to `isolf lipids`.
pub(crate) fn enrich_build_error(force_field: &ForceField, err: BuildError) -> Box<dyn Error> {
    if let BuildError::UnknownLipid(ref name) = err {
        let suggestions = suggest_lipids(force_field, name);
        if suggestions.is_empty() {
            return format!(
                "unknown lipid '{name}' (run 'isolf lipids' to list all {} supported lipids)",
                force_field.lipids.len()
            )
            .into();
        }
        return format!(
            "unknown lipid '{name}' — did you mean {}? (run 'isolf lipids' to list all {})",
            suggestions.join(", "),
            force_field.lipids.len()
        )
        .into();
    }
    err.into()
}

/// Suggest up to 5 closest lipid names within Levenshtein edit distance 2.
fn suggest_lipids<'a>(force_field: &'a ForceField, query: &str) -> Vec<&'a str> {
    let q = query.to_ascii_uppercase();
    let mut ranked: Vec<(usize, &str)> = force_field
        .lipids
        .iter()
        .filter_map(|l| {
            let d = edit_distance(&q, &l.name);
            (d <= 2).then_some((d, l.name.as_str()))
        })
        .collect();
    ranked.sort_by_key(|&(d, name)| (d, name));
    if let Some(&(best, _)) = ranked.first() {
        ranked.retain(|&(d, _)| d == best);
    }
    ranked.into_iter().take(5).map(|(_, n)| n).collect()
}

fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut curr = vec![0; b.len() + 1];
    for (i, &ca) in a.iter().enumerate() {
        curr[0] = i + 1;
        for (j, &cb) in b.iter().enumerate() {
            let cost = usize::from(ca != cb);
            curr[j + 1] = (prev[j + 1] + 1).min(curr[j] + 1).min(prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[b.len()]
}

/// Print the 7×5 iSoLF lipid table directly from the force field.
pub(crate) fn print_lipids(ff: &ForceField, color: bool, unicode: bool) {
    let paint = |s: &str, code: &str| {
        if color {
            format!("\x1b[{code}m{s}\x1b[0m")
        } else {
            s.to_string()
        }
    };
    let tails = [
        ("DL", "lauroyl (12:0)"),
        ("DM", "myristoyl (14:0)"),
        ("DP", "palmitoyl (16:0)"),
        ("DS", "stearoyl (18:0)"),
        ("DO", "dioleoyl (18:1)"),
        ("PO", "16:0 / 18:1"),
        ("SO", "18:0 / 18:1"),
    ];
    let heads = ["PC", "PE", "PG", "PS", "PA"];
    let sep = if unicode { "×" } else { "x" };

    println!(
        "{} ({} species: {} {sep} {}, named <tail><head>)\n",
        paint("iSoLF lipids", "1;36"),
        ff.lipids.len(),
        tails.len(),
        heads.len(),
    );

    let head_cols: Vec<String> = heads
        .iter()
        .map(|h| {
            let sample = format!("PO{h}");
            let q: f64 = ff
                .lipid(&sample)
                .map(|l| {
                    l.beads
                        .iter()
                        .filter_map(|b| ff.bead_type(b))
                        .map(|b| b.charge)
                        .sum()
                })
                .unwrap_or(0.0);
            let q_label = if q.abs() < 0.01 {
                "0".to_string()
            } else {
                format!("{q:+.0}")
            };
            format!("{h} ({q_label})")
        })
        .collect();

    println!(
        "  {}  {}  {}  {}  {}  {}  {}  {}",
        paint(&format!("{:<4}", "Tail"), "1"),
        paint(&format!("{:<18}", "Acyl chains"), "1"),
        paint(&format!("{:<8}", head_cols[0]), "1"),
        paint(&format!("{:<8}", head_cols[1]), "1"),
        paint(&format!("{:<8}", head_cols[2]), "1"),
        paint(&format!("{:<8}", head_cols[3]), "1"),
        paint(&format!("{:<8}", head_cols[4]), "1"),
        paint("Beads", "1"),
    );

    for (tail, desc) in tails {
        let names: Vec<String> = heads.iter().map(|h| format!("{tail}{h}")).collect();
        let pc_beads = ff.lipid(&names[0]).map_or(0, |l| l.beads.len());
        let pa_beads = ff.lipid(&names[4]).map_or(0, |l| l.beads.len());
        println!(
            "  {}  {:<18}  {:<8}  {:<8}  {:<8}  {:<8}  {:<8}  {}",
            paint(&format!("{tail:<4}"), "36"),
            desc,
            names[0],
            names[1],
            names[2],
            names[3],
            names[4],
            paint(&format!("{pc_beads} (PA: {pa_beads})"), "2"),
        );
    }
    println!(
        "\n  {}  isolf membrane -L \"POPC=3,POPS=1\" -n 256 --top -o my_membrane",
        paint("Example:", "2")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn suggests_close_lipid_names() {
        let ff = ForceField::isolf();
        let hits = suggest_lipids(&ff, "POPX");
        assert!(hits.contains(&"POPC"));
        assert!(hits.contains(&"POPS"));
        assert!(suggest_lipids(&ff, "ZZZZZZ").is_empty());
    }
}
