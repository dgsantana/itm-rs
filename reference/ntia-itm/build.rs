//! Compiles NTIA's C++ ITM (`cpp/`, `include/`) into a static library for the reference bindings.

fn main() {
    let mut build = cc::Build::new();
    build.cpp(true).include("include").define("DLLEXPORT", Some(""));
    for entry in std::fs::read_dir("cpp").expect("cpp/") {
        let path = entry.expect("an entry").path();
        if path.extension().is_some_and(|extension| extension == "cpp") {
            println!("cargo:rerun-if-changed={}", path.display());
            build.file(path);
        }
    }
    build.compile("ntia_itm");
}
