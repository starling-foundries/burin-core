//! Freeze burin's own point lookups, forward projections and nuclei bit for bit
//! (`points_burin.json`), and list every disagreement with the reference fixtures
//! (`reference_disagreements.json`). Run `cargo run -p burin-core --release --example freeze_points`;
//! the results are committed and checked exactly by `tests/determinism.rs` and the parity tests.
//! A change to either file is a change to what the crate computes.

#[path = "../tests/support/goldens.rs"]
#[allow(dead_code)]
mod goldens;

use std::path::PathBuf;

fn main() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let own = goldens::self_goldens();
    std::fs::write(dir.join("points_burin.json"), serde_json::to_string(&own).unwrap() + "\n").unwrap();
    let diff = goldens::reference_disagreements(&dir);
    std::fs::write(dir.join("reference_disagreements.json"), serde_json::to_string_pretty(&diff).unwrap() + "\n").unwrap();
    for (name, v) in diff.as_object().unwrap() {
        println!("{name}: {} lookups, {} forward projections, {} nuclei differ from the reference",
                 v["lookups"].as_array().unwrap().len(), v["forward"].as_array().unwrap().len(), v["nuclei"].as_array().unwrap().len());
    }
    let n: usize = own["profiles"].as_object().unwrap().values().map(|p| p["points"].as_array().unwrap().len()).sum();
    println!("wrote points_burin.json ({n} points) and reference_disagreements.json");
}
