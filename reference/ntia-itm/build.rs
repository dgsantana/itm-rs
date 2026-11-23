//! Compiles NTIA's C++ ITM (`cpp/`, `include/`) into a static library for the reference bindings.
//!
//! NTIA's sources are kept exactly as published; they are written for MSVC, so the build compiles copies in `OUT_DIR`
//! with two portable edits: `#include "..\include\x.h"` becomes `#include "x.h"`, and `__declspec(dllexport)` is defined
//! away on compilers that do not know it (`DLLEXPORT` stays `extern "C"`).

use std::{env, fs, path::PathBuf};

fn main() {
    let out = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));
    let sources = out.join("ntia-itm");
    fs::create_dir_all(&sources).expect("a source folder in OUT_DIR");

    for folder in ["cpp", "include"] {
        println!("cargo:rerun-if-changed={folder}");
        for entry in fs::read_dir(folder).unwrap_or_else(|error| panic!("{folder}/: {error}")) {
            let path = entry.expect("an entry").path();
            let text = fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
            let portable = text.replace("#include \"..\\include\\", "#include \"");
            fs::write(sources.join(path.file_name().expect("a file name")), portable).expect("a copy in OUT_DIR");
        }
    }

    let mut build = cc::Build::new();
    build.cpp(true).include(&sources).warnings(false);
    if !build.get_compiler().is_like_msvc() {
        build.define("__declspec(x)", Some(""));
    }
    for entry in fs::read_dir("cpp").expect("cpp/") {
        let name = entry.expect("an entry").file_name();
        if name.to_string_lossy().ends_with(".cpp") {
            build.file(sources.join(name));
        }
    }
    build.compile("ntia_itm");
}
