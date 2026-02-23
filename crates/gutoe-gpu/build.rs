// GUTOE GPU build script — compile CUDA/HIP kernels if the feature is enabled
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

#[cfg(any(feature = "cuda", feature = "rocm"))]
use std::env;
#[cfg(any(feature = "cuda", feature = "rocm"))]
use std::path::PathBuf;
#[cfg(any(feature = "cuda", feature = "rocm"))]
use std::process::Command;

fn main() {
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // ── CUDA (NVIDIA) ──────────────────────────────────────────────────────
    #[cfg(feature = "cuda")]
    {
        let kernel_src = "kernels/cuda_main.cu";
        let obj = out_dir.join("cuda_main.o");
        let lib = out_dir.join("libschrodinger.a");

        // Detect GPU architecture from CUDA_ARCH env, default to sm_89 (RTX 4070 Ti / Ada)
        let arch = env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_89".to_string());

        println!("cargo:warning=Compiling CUDA kernels for {arch}...");

        // Find nvcc: honour NVCC env, then /usr/local/cuda/bin, then PATH
        let nvcc = env::var("NVCC").unwrap_or_else(|_| {
            let default = "/usr/local/cuda/bin/nvcc";
            if std::path::Path::new(default).exists() { default.to_string() }
            else { "nvcc".to_string() }
        });

        let status = Command::new(&nvcc)
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
        println!("cargo:rustc-link-lib=stdc++");   // nvcc C++ runtime symbols

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

        // Find hipcc: honour HIPCC env, then /opt/rocm/bin, then PATH
        let hipcc = env::var("HIPCC").unwrap_or_else(|_| {
            let default = "/opt/rocm/bin/hipcc";
            if std::path::Path::new(default).exists() { default.to_string() }
            else { "hipcc".to_string() }
        });

        // Optional offload arch (e.g. gfx1151 for Strix Halo)
        let mut hipcc_args: Vec<String> = vec![
            "-O3".into(),
            "-fPIC".into(),
            "-D__HIP_PLATFORM_AMD__".into(),
        ];
        if let Ok(gfx) = env::var("GFX_ARCH") {
            hipcc_args.push(format!("--offload-arch={gfx}"));
        }
        hipcc_args.extend(["-c".into(), kernel_src.into(), "-o".into(),
                            obj.to_str().unwrap().to_string()]);

        let status = Command::new(&hipcc)
            .args(&hipcc_args)
            .status()
            .expect("hipcc not found — install ROCm or set HIPCC=/path/to/hipcc");

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
