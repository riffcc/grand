use std::fs;
use std::path::Path;

#[derive(Default, Debug)]
struct Stats {
    rows: usize,
    mad_sum: f64,
    mad_max: f64,
    transfer_abs_sum: f64,
    transfer_abs_max: f64,
    gpu_ms_sum: f64,
    gpu_ms_rows: usize,
    cpu_ms_sum: f64,
    cpu_ms_rows: usize,
}

fn main() {
    let dir = std::env::var("BH_RENDER_DIR").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let root = Path::new(&dir);
    let mut per_case: Vec<(String, String, Stats)> = Vec::new();
    let mut header_only: Vec<String> = Vec::new();

    let entries = match fs::read_dir(root) {
        Ok(e) => e,
        Err(e) => {
            eprintln!("cannot read {}: {e}", root.display());
            std::process::exit(1);
        }
    };

    for ent in entries.flatten() {
        let p = ent.path();
        let Some(name) = p.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        if !name.starts_with("transfer_parity_") || !name.ends_with(".csv") {
            continue;
        }
        let stem = name
            .trim_start_matches("transfer_parity_")
            .trim_end_matches(".csv");
        let (view, backend) = parse_case(stem);
        let mut stats = Stats::default();
        let Ok(text) = fs::read_to_string(&p) else {
            continue;
        };
        let mut lines = text.lines();
        let Some(header) = lines.next() else {
            continue;
        };
        let header_cols: Vec<&str> = header.split(',').collect();
        let idx_mad = header_cols.iter().position(|c| *c == "mad");
        let idx_transfer = header_cols
            .iter()
            .position(|c| *c == "transfer_delta_parity_abs");
        let idx_gpu_ms = header_cols.iter().position(|c| *c == "gpu_ms");
        let idx_cpu_ms = header_cols.iter().position(|c| *c == "cpu_ms");

        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split(',').collect();
            if let Some(i) = idx_mad {
                if let Some(v) = cols.get(i).and_then(|s| s.parse::<f64>().ok()) {
                    stats.mad_sum += v;
                    stats.mad_max = stats.mad_max.max(v);
                }
            }
            if let Some(i) = idx_transfer {
                if let Some(v) = cols.get(i).and_then(|s| s.parse::<f64>().ok()) {
                    stats.transfer_abs_sum += v;
                    stats.transfer_abs_max = stats.transfer_abs_max.max(v);
                }
            }
            if let Some(i) = idx_gpu_ms {
                if let Some(v) = cols.get(i).and_then(|s| s.parse::<f64>().ok()) {
                    stats.gpu_ms_sum += v;
                    stats.gpu_ms_rows += 1;
                }
            }
            if let Some(i) = idx_cpu_ms {
                if let Some(v) = cols.get(i).and_then(|s| s.parse::<f64>().ok()) {
                    stats.cpu_ms_sum += v;
                    stats.cpu_ms_rows += 1;
                }
            }
            stats.rows += 1;
        }
        if stats.rows > 0 {
            per_case.push((view, backend, stats));
        } else {
            header_only.push(name.to_string());
        }
    }

    per_case.sort_by(|a, b| (a.0.as_str(), a.1.as_str()).cmp(&(b.0.as_str(), b.1.as_str())));
    let out = root.join("parity_dashboard.md");
    let mut md = String::new();
    md.push_str("# CUDA/CPU Transfer Parity Dashboard\n\n");
    md.push_str("| view | backend | rows | mean MAD | max MAD | mean |Δtransfer| | max |Δtransfer| | mean gpu_ms | mean cpu_ms | speedup (cpu/gpu) |\n");
    md.push_str("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for (view, backend, s) in &per_case {
        let n = s.rows.max(1) as f64;
        let mad = s.mad_sum / n;
        let td = s.transfer_abs_sum / n;
        let gpu_ms = if s.gpu_ms_rows > 0 {
            s.gpu_ms_sum / s.gpu_ms_rows as f64
        } else {
            f64::NAN
        };
        let cpu_ms = if s.cpu_ms_rows > 0 {
            s.cpu_ms_sum / s.cpu_ms_rows as f64
        } else {
            f64::NAN
        };
        let speedup = if gpu_ms.is_finite() && cpu_ms.is_finite() && gpu_ms > 1e-9 {
            cpu_ms / gpu_ms
        } else {
            f64::NAN
        };
        md.push_str(&format!(
            "| {view} | {backend} | {} | {:.6} | {:.6} | {:.6} | {:.6} | {} | {} | {} |\n",
            s.rows,
            mad,
            s.mad_max,
            td,
            s.transfer_abs_max,
            if gpu_ms.is_finite() {
                format!("{gpu_ms:.3}")
            } else {
                "-".to_string()
            },
            if cpu_ms.is_finite() {
                format!("{cpu_ms:.3}")
            } else {
                "-".to_string()
            },
            if speedup.is_finite() {
                format!("{speedup:.2}x")
            } else {
                "-".to_string()
            }
        ));
    }
    if !header_only.is_empty() {
        md.push_str("\n## Header-Only Parity CSVs\n\n");
        md.push_str("These files contain only headers (no data rows). They usually indicate an aborted run or missing GPU backend at runtime.\n\n");
        for name in &header_only {
            md.push_str(&format!("- `{name}`\n"));
        }
    }
    fs::write(&out, md).expect("write dashboard");
    println!("wrote {}", out.display());
    if per_case.is_empty() {
        eprintln!("no parity data rows found; dashboard includes header-only files for diagnosis");
        std::process::exit(2);
    }
}

fn parse_case(stem: &str) -> (String, String) {
    for tag in ["cuda", "rocm", "hip", "cpu"] {
        let suffix = format!("_{tag}");
        if let Some(view) = stem.strip_suffix(&suffix) {
            return (view.to_string(), tag.to_string());
        }
    }
    (stem.to_string(), "unknown".to_string())
}
