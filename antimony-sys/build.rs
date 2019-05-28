extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-link-lib=antimony");

    let path = PathBuf::from(env::var("OUT_DIR").unwrap()).join("bindings.rs");

    bindgen::Builder::default()
        .header("wrapper.h")
        .default_enum_style(bindgen::EnumVariation::Rust)
        .blacklist_type("var_type")
        .blacklist_type("const_type")
        .blacklist_type("deletion_type")
        .blacklist_type("distribution_type")
        .blacklist_function("freeAll")
        .generate_comments(true)
        .generate()
        .expect("Failed to generate bindings")
        .write_to_file(path)
        .expect("Failed to write bindings");
}
