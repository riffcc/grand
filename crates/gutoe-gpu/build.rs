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
        // Detect GPU architecture from CUDA_ARCH env, default to sm_89 (RTX 4070 Ti / Ada)
        let arch = env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_89".to_string());

        println!("cargo:warning=Compiling CUDA kernels for {arch}...");

        // Find nvcc: honour NVCC env, then /usr/local/cuda/bin, then PATH
        let nvcc = env::var("NVCC").unwrap_or_else(|_| {
            let default = "/usr/local/cuda/bin/nvcc";
            if std::path::Path::new(default).exists() { default.to_string() }
            else { "nvcc".to_string() }
        });

        // Compile all kernel sources to object files
        let sources = [
            ("kernels/cuda_main.cu", out_dir.join("cuda_main.o")),
            ("kernels/tracer.cu",    out_dir.join("tracer_gpu.o")),
        ];
        let mut obj_paths: Vec<String> = Vec::new();

        for (src, obj) in &sources {
            let status = Command::new(&nvcc)
                .args(&[
                    "-O3",
                    &format!("-arch={arch}"),
                    "-Xcompiler", "-fPIC",
                    "-c", src,
                    "-o", obj.to_str().unwrap(),
                ])
                .status()
                .unwrap_or_else(|_| panic!("nvcc not found — install CUDA toolkit or unset --features cuda"));

            if !status.success() {
                panic!("nvcc compilation failed for {src}");
            }
            obj_paths.push(obj.to_str().unwrap().to_string());
        }

        // Archive all objects into one static lib
        let lib = out_dir.join("libschrodinger.a");
        let mut ar_args = vec!["rcs".to_string(), lib.to_str().unwrap().to_string()];
        ar_args.extend(obj_paths);
        Command::new("ar")
            .args(&ar_args)
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
        println!("cargo:rerun-if-changed=kernels/tracer.cu");
    }

    // ── ROCm / HIP (AMD) ───────────────────────────────────────────────────
    #[cfg(feature = "rocm")]
    {
        println!("cargo:warning=Compiling HIP kernels for ROCm...");

        // Find hipcc: honour HIPCC env, then /opt/rocm/bin, then PATH
        let hipcc = env::var("HIPCC").unwrap_or_else(|_| {
            let default = "/opt/rocm/bin/hipcc";
            if std::path::Path::new(default).exists() { default.to_string() }
            else { "hipcc".to_string() }
        });

        // Optional offload arch (e.g. gfx1151 for Strix Halo)
        let mut base_args: Vec<String> = vec![
            "-O3".into(),
            "-fPIC".into(),
            "-D__HIP_PLATFORM_AMD__".into(),
        ];
        // Some HIP header/toolchain mixes need this macro during host parsing,
        // while others reject the wrong value for a given target arch. We try
        // a small candidate set and keep the first one that compiles.
        let gfx_arch = env::var("GFX_ARCH").ok();
        // Prefer ROCm include tree when available to avoid mixed system headers.
        let rocm_inc = env::var("ROCM_PATH")
            .map(|r| format!("{r}/include"))
            .unwrap_or_else(|_| "/opt/rocm/include".to_string());
        base_args.push(format!("-I{rocm_inc}"));
        if let Ok(gfx) = env::var("GFX_ARCH") {
            base_args.push(format!("--offload-arch={gfx}"));
        }

        // Compile all kernel sources (same .cu sources work with HIP)
        let sources = [
            ("kernels/cuda_main.cu", out_dir.join("hip_main.o")),
            ("kernels/tracer.cu",    out_dir.join("hip_tracer.o")),
        ];
        let mut selected_wavefront: Option<String> = None;
        let wavefront_candidates: Vec<Option<String>> = if let Ok(v) = env::var("HIP_WAVEFRONT_SIZE") {
            vec![Some(v)]
        } else if gfx_arch
            .as_deref()
            .is_some_and(|g| g.starts_with("gfx11"))
        {
            vec![Some("32".to_string()), None, Some("64".to_string())]
        } else {
            vec![None, Some("64".to_string()), Some("32".to_string())]
        };

        let mut obj_paths: Vec<String> = Vec::new();
        let mut compiled = false;
        for wf in wavefront_candidates {
            obj_paths.clear();
            let mut ok_all = true;
            for (src, obj) in &sources {
                let mut args = base_args.clone();
                if let Some(ref wf_val) = wf {
                    args.push(format!("-D__AMDGCN_WAVEFRONT_SIZE={wf_val}"));
                }
                args.extend([
                    "-c".into(),
                    src.to_string(),
                    "-o".into(),
                    obj.to_str().unwrap().to_string(),
                ]);
                let status = Command::new(&hipcc)
                    .args(&args)
                    .status()
                    .expect("hipcc not found — install ROCm or set HIPCC=/path/to/hipcc");
                if !status.success() {
                    ok_all = false;
                    break;
                }
                obj_paths.push(obj.to_str().unwrap().to_string());
            }
            if ok_all {
                compiled = true;
                selected_wavefront = wf;
                break;
            }
        }
        if !compiled {
            panic!("hipcc compilation failed for kernels/cuda_main.cu");
        }
        println!(
            "cargo:warning=HIP wavefront={} (GFX_ARCH={})",
            selected_wavefront.as_deref().unwrap_or("auto"),
            gfx_arch.as_deref().unwrap_or("auto")
        );

        let lib = out_dir.join("libschrodinger.a");
        let mut ar_args = vec!["rcs".to_string(), lib.to_str().unwrap().to_string()];
        ar_args.extend(obj_paths);
        Command::new("ar").args(&ar_args).status().expect("ar failed");

        println!("cargo:rustc-link-search=native={}", out_dir.display());
        println!("cargo:rustc-link-lib=static=schrodinger");

        let rocm_lib = env::var("ROCM_PATH")
            .map(|r| format!("{r}/lib"))
            .unwrap_or_else(|_| "/opt/rocm/lib".to_string());
        println!("cargo:rustc-link-search=native={rocm_lib}");
        println!("cargo:rustc-link-lib=amdhip64");
        println!("cargo:rerun-if-changed=kernels/cuda_main.cu");
        println!("cargo:rerun-if-changed=kernels/tracer.cu");
    }
}
