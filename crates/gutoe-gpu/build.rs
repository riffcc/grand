// GUTOE GPU build script — compile CUDA/HIP kernels if the feature is enabled
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // ── CUDA (NVIDIA) ──────────────────────────────────────────────────────
    #[cfg(feature = "cuda")]
    {
        let kernel_src = "kernels/cuda_main.cu";
        let obj = out_dir.join("cuda_main.o");
        let lib = out_dir.join("libschrodinger.a");

        // Detect GPU architecture from CUDA_ARCH env, default to sm_86 (3070 Ti)
        let arch = env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_86".to_string());

        println!("cargo:warning=Compiling CUDA kernels for {arch}...");

        let status = Command::new("nvcc")
            .args(&[
                "-O3",
                &format!("-arch={arch}"),
                "-Xcompiler", "-fPIC",
                "-c", kernel_src,
                "-o", obj.to_str().unwrap(),
            ])
            .status()
            .expect("nvcc not found — install CUDA toolkit or unset --features cuda");

        if !status.success() {
            panic!("nvcc compilation failed");
        }

        Command::new("ar")
            .args(&["rcs", lib.to_str().unwrap(), obj.to_str().unwrap()])
            .status()
            .expect("ar failed");

        println!("cargo:rustc-link-search=native={}", out_dir.display());
        println!("cargo:rustc-link-lib=static=schrodinger");
        println!("cargo:rustc-link-lib=cudart");

        // Find CUDA lib path
        let cuda_lib = env::var("CUDA_ROOT")
            .map(|r| format!("{r}/lib64"))
            .unwrap_or_else(|_| "/usr/local/cuda/lib64".to_string());
        println!("cargo:rustc-link-search=native={cuda_lib}");
        println!("cargo:rerun-if-changed=kernels/cuda_main.cu");
    }

    // ── ROCm / HIP (AMD) ───────────────────────────────────────────────────
    #[cfg(feature = "rocm")]
    {
        let kernel_src = "kernels/cuda_main.cu"; // same source, HIP is compatible
        let obj = out_dir.join("hip_main.o");
        let lib = out_dir.join("libschrodinger.a");

        println!("cargo:warning=Compiling HIP kernels for ROCm...");

        let status = Command::new("hipcc")
            .args(&[
                "-O3",
                "-Xcompiler", "-fPIC",
                "-c", kernel_src,
                "-o", obj.to_str().unwrap(),
            ])
            .status()
            .expect("hipcc not found — install ROCm or unset --features rocm");

        if !status.success() {
            panic!("hipcc compilation failed");
        }

        Command::new("ar")
            .args(&["rcs", lib.to_str().unwrap(), obj.to_str().unwrap()])
            .status()
            .expect("ar failed");

        println!("cargo:rustc-link-search=native={}", out_dir.display());
        println!("cargo:rustc-link-lib=static=schrodinger");

        let rocm_lib = env::var("ROCM_PATH")
            .map(|r| format!("{r}/lib"))
            .unwrap_or_else(|_| "/opt/rocm/lib".to_string());
        println!("cargo:rustc-link-search=native={rocm_lib}");
        println!("cargo:rustc-link-lib=amdhip64");
        println!("cargo:rerun-if-changed=kernels/cuda_main.cu");
    }
}
