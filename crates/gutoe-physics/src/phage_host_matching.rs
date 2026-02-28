/*!
 * Phage-host matching lane (reduced-order binding/lysis proxy).
 *
 * Purpose:
 * - Rank phage candidates against resistant bacterial strains by
 *   receptor-binding energetics plus takeover efficiency.
 * - Explicitly model resistance-bypass behavior: beta-lactamase class does not
 *   enter the phage binding path.
 *
 * This is simulation infrastructure, not clinical guidance.
 */

use crate::chemical_thermo::{AVOGADRO, R_GAS_J_MOL_K};
use crate::{ALPHA_LEADING_ORDER, C, HBAR};
use serde_json::Value;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReceptorKind {
    LamB,
    OmpK35,
    OmpK36,
    FhuA,
    LpsCore,
    TypeIvPilus,
}

#[derive(Clone, Copy, Debug)]
pub struct ReceptorProfile {
    pub lamb: f64,
    pub ompk35: f64,
    pub ompk36: f64,
    pub fhua: f64,
    pub lps_core: f64,
    pub type_iv_pilus: f64,
}

#[derive(Clone, Debug)]
pub struct BacterialStrainSpec {
    pub name: String,
    pub species: String,
    pub resistance_marker: String,
    pub receptor_profile: ReceptorProfile,
}

#[derive(Clone, Debug)]
pub struct PhageSpec {
    pub name: String,
    pub family: String,
    pub primary_receptor: ReceptorKind,
    pub secondary_receptor: Option<ReceptorKind>,
    pub secondary_weight: f64,
    pub ionic_contact_count: f64,
    pub hbond_contact_count: f64,
    pub hydrophobic_area_a2: f64,
    pub conformational_entropy_penalty: f64,
    pub host_takeover_efficiency: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct PhageMatchingCoefficients {
    pub ionic_distance_nm: f64,
    pub hbond_distance_nm: f64,
    pub active_site_dielectric: f64,
    pub hbond_charge_product: f64,
    pub hydrophobic_coeff_kj_per_a2: f64,
    pub mismatch_penalty_kj: f64,
    pub baseline_offset_kj: f64,
    pub local_effective_concentration_nanomolar: f64,
}

impl Default for PhageMatchingCoefficients {
    fn default() -> Self {
        Self {
            ionic_distance_nm: 0.36,
            hbond_distance_nm: 0.31,
            active_site_dielectric: 30.0,
            hbond_charge_product: 0.18,
            hydrophobic_coeff_kj_per_a2: 0.0065,
            mismatch_penalty_kj: 8.5,
            baseline_offset_kj: -15.0,
            local_effective_concentration_nanomolar: 100.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct PhageHostPairResult {
    pub strain_name: String,
    pub strain_species: String,
    pub resistance_marker: String,
    pub phage_name: String,
    pub phage_family: String,
    pub receptor_match_score: f64,
    pub qed_ionic_floor_kj_mol: f64,
    pub qed_hbond_floor_kj_mol: f64,
    pub qed_floor_total_kj_mol: f64,
    pub residual_modeled_total_kj_mol: f64,
    pub predicted_delta_g_kj_mol: f64,
    pub predicted_kd_nanomolar: f64,
    pub attachment_prob: f64,
    pub lysis_potential_score: f64,
}

#[derive(Clone, Debug)]
pub struct StrainBestMatch {
    pub strain_name: String,
    pub resistance_marker: String,
    pub best_phage_name: String,
    pub best_lysis_score: f64,
    pub best_predicted_kd_nanomolar: f64,
}

#[derive(Clone, Debug)]
pub struct PhageMatchingPanel {
    pub rows: Vec<PhageHostPairResult>,
    pub best_by_strain: Vec<StrainBestMatch>,
    pub resistance_independence_probe_abs_delta: f64,
    pub mean_best_lysis_score: f64,
}

fn receptor_expression(profile: ReceptorProfile, receptor: ReceptorKind) -> f64 {
    match receptor {
        ReceptorKind::LamB => profile.lamb,
        ReceptorKind::OmpK35 => profile.ompk35,
        ReceptorKind::OmpK36 => profile.ompk36,
        ReceptorKind::FhuA => profile.fhua,
        ReceptorKind::LpsCore => profile.lps_core,
        ReceptorKind::TypeIvPilus => profile.type_iv_pilus,
    }
    .clamp(0.0, 1.0)
}

fn qed_contact_energy_kj_mol(charge_product: f64, distance_nm: f64, dielectric: f64) -> f64 {
    let q = charge_product.abs();
    let r_m = distance_nm.max(1.0e-6) * 1.0e-9;
    let eps = dielectric.max(1.0);
    let per_molecule_j = -(q * ALPHA_LEADING_ORDER * HBAR * C) / (eps * r_m);
    per_molecule_j * AVOGADRO / 1000.0
}

fn kd_nanomolar_from_delta_g(delta_g_kj_mol: f64, temperature_k: f64) -> f64 {
    let exponent = delta_g_kj_mol * 1000.0 / (R_GAS_J_MOL_K * temperature_k.max(1.0));
    exponent.exp() * 1.0e9
}

fn pairwise_score(
    strain: &BacterialStrainSpec,
    phage: &PhageSpec,
    temperature_k: f64,
    c: PhageMatchingCoefficients,
) -> PhageHostPairResult {
    let primary = receptor_expression(strain.receptor_profile, phage.primary_receptor);
    let secondary = phage
        .secondary_receptor
        .map(|r| receptor_expression(strain.receptor_profile, r))
        .unwrap_or(0.0);
    let match_score = (primary + phage.secondary_weight.max(0.0) * secondary).clamp(0.0, 1.0);

    let ionic_floor = phage.ionic_contact_count.max(0.0)
        * match_score
        * qed_contact_energy_kj_mol(1.0, c.ionic_distance_nm, c.active_site_dielectric);
    let hbond_floor = phage.hbond_contact_count.max(0.0)
        * match_score
        * qed_contact_energy_kj_mol(
            c.hbond_charge_product,
            c.hbond_distance_nm,
            c.active_site_dielectric + 2.0,
        );
    let qed_total = ionic_floor + hbond_floor;

    let residual = -c.hydrophobic_coeff_kj_per_a2
        * phage.hydrophobic_area_a2.max(0.0)
        * match_score
        + c.mismatch_penalty_kj * (1.0 - match_score)
        + phage.conformational_entropy_penalty.max(0.0)
        + c.baseline_offset_kj;

    let predicted_delta_g = qed_total + residual;
    let kd_nanomolar = kd_nanomolar_from_delta_g(predicted_delta_g, temperature_k).max(1.0e-9);
    let attachment_prob = c.local_effective_concentration_nanomolar
        / (c.local_effective_concentration_nanomolar + kd_nanomolar);
    let lysis_score = (attachment_prob * phage.host_takeover_efficiency.max(0.0)).clamp(0.0, 1.0);

    PhageHostPairResult {
        strain_name: strain.name.clone(),
        strain_species: strain.species.clone(),
        resistance_marker: strain.resistance_marker.clone(),
        phage_name: phage.name.clone(),
        phage_family: phage.family.clone(),
        receptor_match_score: match_score,
        qed_ionic_floor_kj_mol: ionic_floor,
        qed_hbond_floor_kj_mol: hbond_floor,
        qed_floor_total_kj_mol: qed_total,
        residual_modeled_total_kj_mol: residual,
        predicted_delta_g_kj_mol: predicted_delta_g,
        predicted_kd_nanomolar: kd_nanomolar,
        attachment_prob,
        lysis_potential_score: lysis_score,
    }
}

fn best_for_strain(rows: &[PhageHostPairResult], strain_name: &str) -> Option<StrainBestMatch> {
    let best = rows
        .iter()
        .filter(|r| r.strain_name == strain_name)
        .max_by(|a, b| {
            a.lysis_potential_score
                .partial_cmp(&b.lysis_potential_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })?
        .clone();
    Some(StrainBestMatch {
        strain_name: best.strain_name.clone(),
        resistance_marker: best.resistance_marker.clone(),
        best_phage_name: best.phage_name.clone(),
        best_lysis_score: best.lysis_potential_score,
        best_predicted_kd_nanomolar: best.predicted_kd_nanomolar,
    })
}

pub fn parse_receptor_kind(s: &str) -> Result<ReceptorKind, String> {
    let key = s.trim().to_ascii_lowercase().replace(['-', '_', ' '], "");
    match key.as_str() {
        "lamb" => Ok(ReceptorKind::LamB),
        "ompk35" => Ok(ReceptorKind::OmpK35),
        "ompk36" => Ok(ReceptorKind::OmpK36),
        "fhua" => Ok(ReceptorKind::FhuA),
        "lpscore" => Ok(ReceptorKind::LpsCore),
        "typeivpilus" => Ok(ReceptorKind::TypeIvPilus),
        _ => Err(format!("unknown receptor kind: {s}")),
    }
}

fn normalize_header(s: &str) -> String {
    s.trim()
        .to_ascii_lowercase()
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .collect::<String>()
}

fn parse_csv_line(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        match ch {
            '"' => {
                if in_quotes && chars.peek() == Some(&'"') {
                    cur.push('"');
                    let _ = chars.next();
                } else {
                    in_quotes = !in_quotes;
                }
            }
            ',' if !in_quotes => {
                out.push(cur.trim().to_string());
                cur.clear();
            }
            _ => cur.push(ch),
        }
    }

    out.push(cur.trim().to_string());
    out
}

fn find_col(headers: &[String], aliases: &[&str]) -> Option<usize> {
    let aliases_norm = aliases
        .iter()
        .map(|a| normalize_header(a))
        .collect::<Vec<_>>();
    headers.iter().position(|h| aliases_norm.iter().any(|a| a == h))
}

fn get_required_cell<'a>(row: &'a [String], idx: usize, label: &str) -> Result<&'a str, String> {
    row.get(idx)
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .ok_or_else(|| format!("missing required field: {label}"))
}

fn get_f64_cell(row: &[String], idx: usize, label: &str) -> Result<f64, String> {
    let s = get_required_cell(row, idx, label)?;
    s.parse::<f64>()
        .map_err(|e| format!("invalid float for {label}: {s} ({e})"))
}

fn csv_rows(data: &str) -> Result<Vec<Vec<String>>, String> {
    let rows = data
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(parse_csv_line)
        .collect::<Vec<_>>();
    if rows.is_empty() {
        return Err("csv has no rows".to_string());
    }
    Ok(rows)
}

fn v_required_f64(obj: &Value, key: &str) -> Result<f64, String> {
    obj.get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("missing or invalid f64 field: {key}"))
}

fn v_required_str(obj: &Value, key: &str) -> Result<String, String> {
    obj.get(key)
        .and_then(Value::as_str)
        .map(|s| s.to_string())
        .ok_or_else(|| format!("missing or invalid string field: {key}"))
}

fn json_array_or_key<'a>(v: &'a Value, key: &str) -> Result<&'a Vec<Value>, String> {
    if let Some(arr) = v.as_array() {
        return Ok(arr);
    }
    v.get(key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("expected JSON array or object with key '{key}'"))
}

pub fn load_strains_from_json_str(data: &str) -> Result<Vec<BacterialStrainSpec>, String> {
    let root: Value = serde_json::from_str(data).map_err(|e| format!("json parse error: {e}"))?;
    let arr = json_array_or_key(&root, "strains")?;

    let mut out = Vec::with_capacity(arr.len());
    for (i, row) in arr.iter().enumerate() {
        let rec = if row.get("receptors").is_some() {
            row.get("receptors").unwrap_or(row)
        } else {
            row
        };
        let strain = BacterialStrainSpec {
            name: v_required_str(row, "name").map_err(|e| format!("row {i}: {e}"))?,
            species: v_required_str(row, "species").map_err(|e| format!("row {i}: {e}"))?,
            resistance_marker: v_required_str(row, "resistance_marker")
                .map_err(|e| format!("row {i}: {e}"))?,
            receptor_profile: ReceptorProfile {
                lamb: v_required_f64(rec, "lamb").map_err(|e| format!("row {i}: {e}"))?,
                ompk35: v_required_f64(rec, "ompk35").map_err(|e| format!("row {i}: {e}"))?,
                ompk36: v_required_f64(rec, "ompk36").map_err(|e| format!("row {i}: {e}"))?,
                fhua: v_required_f64(rec, "fhua").map_err(|e| format!("row {i}: {e}"))?,
                lps_core: v_required_f64(rec, "lps_core")
                    .or_else(|_| v_required_f64(rec, "lpscore"))
                    .map_err(|e| format!("row {i}: {e}"))?,
                type_iv_pilus: v_required_f64(rec, "type_iv_pilus")
                    .or_else(|_| v_required_f64(rec, "typeivpilus"))
                    .map_err(|e| format!("row {i}: {e}"))?,
            },
        };
        out.push(strain);
    }
    Ok(out)
}

pub fn load_phages_from_json_str(data: &str) -> Result<Vec<PhageSpec>, String> {
    let root: Value = serde_json::from_str(data).map_err(|e| format!("json parse error: {e}"))?;
    let arr = json_array_or_key(&root, "phages")?;

    let mut out = Vec::with_capacity(arr.len());
    for (i, row) in arr.iter().enumerate() {
        let primary = v_required_str(row, "primary_receptor").map_err(|e| format!("row {i}: {e}"))?;
        let secondary = row
            .get("secondary_receptor")
            .and_then(Value::as_str)
            .map(parse_receptor_kind)
            .transpose()
            .map_err(|e| format!("row {i}: {e}"))?;

        let phage = PhageSpec {
            name: v_required_str(row, "name").map_err(|e| format!("row {i}: {e}"))?,
            family: v_required_str(row, "family").map_err(|e| format!("row {i}: {e}"))?,
            primary_receptor: parse_receptor_kind(&primary)
                .map_err(|e| format!("row {i}: {e}"))?,
            secondary_receptor: secondary,
            secondary_weight: v_required_f64(row, "secondary_weight")
                .map_err(|e| format!("row {i}: {e}"))?,
            ionic_contact_count: v_required_f64(row, "ionic_contact_count")
                .map_err(|e| format!("row {i}: {e}"))?,
            hbond_contact_count: v_required_f64(row, "hbond_contact_count")
                .map_err(|e| format!("row {i}: {e}"))?,
            hydrophobic_area_a2: v_required_f64(row, "hydrophobic_area_a2")
                .map_err(|e| format!("row {i}: {e}"))?,
            conformational_entropy_penalty: v_required_f64(row, "conformational_entropy_penalty")
                .map_err(|e| format!("row {i}: {e}"))?,
            host_takeover_efficiency: v_required_f64(row, "host_takeover_efficiency")
                .map_err(|e| format!("row {i}: {e}"))?,
        };
        out.push(phage);
    }
    Ok(out)
}

pub fn load_strains_from_csv_str(data: &str) -> Result<Vec<BacterialStrainSpec>, String> {
    let rows = csv_rows(data)?;
    let headers = rows[0].iter().map(|h| normalize_header(h)).collect::<Vec<_>>();

    let idx_name = find_col(&headers, &["name"])
        .ok_or_else(|| "missing 'name' column".to_string())?;
    let idx_species = find_col(&headers, &["species"])
        .ok_or_else(|| "missing 'species' column".to_string())?;
    let idx_res = find_col(&headers, &["resistance_marker", "resistance"])
        .ok_or_else(|| "missing 'resistance_marker' column".to_string())?;
    let idx_lamb = find_col(&headers, &["lamb"])
        .ok_or_else(|| "missing 'lamb' column".to_string())?;
    let idx_ompk35 = find_col(&headers, &["ompk35"])
        .ok_or_else(|| "missing 'ompk35' column".to_string())?;
    let idx_ompk36 = find_col(&headers, &["ompk36"])
        .ok_or_else(|| "missing 'ompk36' column".to_string())?;
    let idx_fhua = find_col(&headers, &["fhua"])
        .ok_or_else(|| "missing 'fhua' column".to_string())?;
    let idx_lps = find_col(&headers, &["lps_core", "lpscore"])
        .ok_or_else(|| "missing 'lps_core' column".to_string())?;
    let idx_pilus = find_col(&headers, &["type_iv_pilus", "typeivpilus"])
        .ok_or_else(|| "missing 'type_iv_pilus' column".to_string())?;

    let mut out = Vec::new();
    for (line_idx, row) in rows.iter().enumerate().skip(1) {
        let mk_err = |e: String| format!("line {}: {}", line_idx + 1, e);
        out.push(BacterialStrainSpec {
            name: get_required_cell(row, idx_name, "name").map_err(mk_err.clone())?.to_string(),
            species: get_required_cell(row, idx_species, "species")
                .map_err(mk_err.clone())?
                .to_string(),
            resistance_marker: get_required_cell(row, idx_res, "resistance_marker")
                .map_err(mk_err.clone())?
                .to_string(),
            receptor_profile: ReceptorProfile {
                lamb: get_f64_cell(row, idx_lamb, "lamb").map_err(mk_err.clone())?,
                ompk35: get_f64_cell(row, idx_ompk35, "ompk35").map_err(mk_err.clone())?,
                ompk36: get_f64_cell(row, idx_ompk36, "ompk36").map_err(mk_err.clone())?,
                fhua: get_f64_cell(row, idx_fhua, "fhua").map_err(mk_err.clone())?,
                lps_core: get_f64_cell(row, idx_lps, "lps_core").map_err(mk_err.clone())?,
                type_iv_pilus: get_f64_cell(row, idx_pilus, "type_iv_pilus").map_err(mk_err)?,
            },
        });
    }
    Ok(out)
}

pub fn load_phages_from_csv_str(data: &str) -> Result<Vec<PhageSpec>, String> {
    let rows = csv_rows(data)?;
    let headers = rows[0].iter().map(|h| normalize_header(h)).collect::<Vec<_>>();

    let idx_name = find_col(&headers, &["name"])
        .ok_or_else(|| "missing 'name' column".to_string())?;
    let idx_family = find_col(&headers, &["family"])
        .ok_or_else(|| "missing 'family' column".to_string())?;
    let idx_pr = find_col(&headers, &["primary_receptor"])
        .ok_or_else(|| "missing 'primary_receptor' column".to_string())?;
    let idx_sr = find_col(&headers, &["secondary_receptor"])
        .ok_or_else(|| "missing 'secondary_receptor' column".to_string())?;
    let idx_sw = find_col(&headers, &["secondary_weight"])
        .ok_or_else(|| "missing 'secondary_weight' column".to_string())?;
    let idx_ic = find_col(&headers, &["ionic_contact_count"])
        .ok_or_else(|| "missing 'ionic_contact_count' column".to_string())?;
    let idx_hc = find_col(&headers, &["hbond_contact_count"])
        .ok_or_else(|| "missing 'hbond_contact_count' column".to_string())?;
    let idx_ha = find_col(&headers, &["hydrophobic_area_a2"])
        .ok_or_else(|| "missing 'hydrophobic_area_a2' column".to_string())?;
    let idx_ce = find_col(&headers, &["conformational_entropy_penalty"])
        .ok_or_else(|| "missing 'conformational_entropy_penalty' column".to_string())?;
    let idx_ht = find_col(&headers, &["host_takeover_efficiency"])
        .ok_or_else(|| "missing 'host_takeover_efficiency' column".to_string())?;

    let mut out = Vec::new();
    for (line_idx, row) in rows.iter().enumerate().skip(1) {
        let mk_err = |e: String| format!("line {}: {}", line_idx + 1, e);
        let primary_raw = get_required_cell(row, idx_pr, "primary_receptor").map_err(mk_err.clone())?;
        let secondary_raw = row
            .get(idx_sr)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .unwrap_or("");

        out.push(PhageSpec {
            name: get_required_cell(row, idx_name, "name").map_err(mk_err.clone())?.to_string(),
            family: get_required_cell(row, idx_family, "family").map_err(mk_err.clone())?.to_string(),
            primary_receptor: parse_receptor_kind(primary_raw).map_err(mk_err.clone())?,
            secondary_receptor: if secondary_raw.is_empty() {
                None
            } else {
                Some(parse_receptor_kind(secondary_raw).map_err(mk_err.clone())?)
            },
            secondary_weight: get_f64_cell(row, idx_sw, "secondary_weight").map_err(mk_err.clone())?,
            ionic_contact_count: get_f64_cell(row, idx_ic, "ionic_contact_count").map_err(mk_err.clone())?,
            hbond_contact_count: get_f64_cell(row, idx_hc, "hbond_contact_count").map_err(mk_err.clone())?,
            hydrophobic_area_a2: get_f64_cell(row, idx_ha, "hydrophobic_area_a2").map_err(mk_err.clone())?,
            conformational_entropy_penalty: get_f64_cell(
                row,
                idx_ce,
                "conformational_entropy_penalty",
            )
            .map_err(mk_err.clone())?,
            host_takeover_efficiency: get_f64_cell(row, idx_ht, "host_takeover_efficiency")
                .map_err(mk_err)?,
        });
    }
    Ok(out)
}

pub fn load_strains_from_path(path: &str) -> Result<Vec<BacterialStrainSpec>, String> {
    let data = std::fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".json") {
        load_strains_from_json_str(&data)
    } else if lower.ends_with(".csv") {
        load_strains_from_csv_str(&data)
    } else {
        load_strains_from_json_str(&data).or_else(|_| load_strains_from_csv_str(&data))
    }
}

pub fn load_phages_from_path(path: &str) -> Result<Vec<PhageSpec>, String> {
    let data = std::fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".json") {
        load_phages_from_json_str(&data)
    } else if lower.ends_with(".csv") {
        load_phages_from_csv_str(&data)
    } else {
        load_phages_from_json_str(&data).or_else(|_| load_phages_from_csv_str(&data))
    }
}

pub fn default_phage_specs() -> Vec<PhageSpec> {
    vec![
        PhageSpec {
            name: "phi_kp_omp".to_string(),
            family: "myoviridae_like".to_string(),
            primary_receptor: ReceptorKind::OmpK36,
            secondary_receptor: Some(ReceptorKind::OmpK35),
            secondary_weight: 0.60,
            ionic_contact_count: 2.8,
            hbond_contact_count: 4.1,
            hydrophobic_area_a2: 360.0,
            conformational_entropy_penalty: 2.6,
            host_takeover_efficiency: 0.88,
        },
        PhageSpec {
            name: "phi_lambda_lamb".to_string(),
            family: "siphoviridae_like".to_string(),
            primary_receptor: ReceptorKind::LamB,
            secondary_receptor: Some(ReceptorKind::LpsCore),
            secondary_weight: 0.35,
            ionic_contact_count: 2.4,
            hbond_contact_count: 3.8,
            hydrophobic_area_a2: 330.0,
            conformational_entropy_penalty: 2.4,
            host_takeover_efficiency: 0.82,
        },
        PhageSpec {
            name: "phi_fhua_spike".to_string(),
            family: "podoviridae_like".to_string(),
            primary_receptor: ReceptorKind::FhuA,
            secondary_receptor: Some(ReceptorKind::LpsCore),
            secondary_weight: 0.30,
            ionic_contact_count: 2.1,
            hbond_contact_count: 3.5,
            hydrophobic_area_a2: 315.0,
            conformational_entropy_penalty: 2.3,
            host_takeover_efficiency: 0.76,
        },
        PhageSpec {
            name: "phi_lps_broad".to_string(),
            family: "myoviridae_like".to_string(),
            primary_receptor: ReceptorKind::LpsCore,
            secondary_receptor: None,
            secondary_weight: 0.0,
            ionic_contact_count: 1.9,
            hbond_contact_count: 3.1,
            hydrophobic_area_a2: 280.0,
            conformational_entropy_penalty: 2.1,
            host_takeover_efficiency: 0.70,
        },
    ]
}

pub fn default_strain_specs() -> Vec<BacterialStrainSpec> {
    vec![
        BacterialStrainSpec {
            name: "kp_ndm1_clinical".to_string(),
            species: "Klebsiella pneumoniae".to_string(),
            resistance_marker: "NDM-1".to_string(),
            receptor_profile: ReceptorProfile {
                lamb: 0.05,
                ompk35: 0.72,
                ompk36: 0.86,
                fhua: 0.32,
                lps_core: 0.90,
                type_iv_pilus: 0.10,
            },
        },
        BacterialStrainSpec {
            name: "kp_kpc_clinical".to_string(),
            species: "Klebsiella pneumoniae".to_string(),
            resistance_marker: "KPC".to_string(),
            receptor_profile: ReceptorProfile {
                lamb: 0.04,
                ompk35: 0.58,
                ompk36: 0.83,
                fhua: 0.30,
                lps_core: 0.86,
                type_iv_pilus: 0.10,
            },
        },
        BacterialStrainSpec {
            name: "ec_tem1_clinical".to_string(),
            species: "Escherichia coli".to_string(),
            resistance_marker: "TEM-1".to_string(),
            receptor_profile: ReceptorProfile {
                lamb: 0.88,
                ompk35: 0.18,
                ompk36: 0.20,
                fhua: 0.73,
                lps_core: 0.82,
                type_iv_pilus: 0.12,
            },
        },
        BacterialStrainSpec {
            name: "pa_mdr_clinical".to_string(),
            species: "Pseudomonas aeruginosa".to_string(),
            resistance_marker: "VIM-2".to_string(),
            receptor_profile: ReceptorProfile {
                lamb: 0.02,
                ompk35: 0.10,
                ompk36: 0.12,
                fhua: 0.18,
                lps_core: 0.64,
                type_iv_pilus: 0.86,
            },
        },
    ]
}

pub fn evaluate_phage_host_matching_panel(
    strains: &[BacterialStrainSpec],
    phages: &[PhageSpec],
    temperature_k: f64,
    coeffs: PhageMatchingCoefficients,
) -> PhageMatchingPanel {
    let mut rows = Vec::new();
    for strain in strains {
        for phage in phages {
            rows.push(pairwise_score(strain, phage, temperature_k, coeffs));
        }
    }

    let best_by_strain = strains
        .iter()
        .filter_map(|s| best_for_strain(&rows, &s.name))
        .collect::<Vec<_>>();

    let mean_best_lysis_score = best_by_strain
        .iter()
        .map(|b| b.best_lysis_score)
        .sum::<f64>()
        / best_by_strain.len().max(1) as f64;

    let resistance_independence_probe_abs_delta = if phages.is_empty() {
        0.0
    } else {
        // Independence probe: identical receptor profile, resistance marker flip.
        let probe_profile = ReceptorProfile {
            lamb: 0.10,
            ompk35: 0.70,
            ompk36: 0.82,
            fhua: 0.30,
            lps_core: 0.88,
            type_iv_pilus: 0.10,
        };
        let probe_a = BacterialStrainSpec {
            name: "probe_a".to_string(),
            species: "Klebsiella pneumoniae".to_string(),
            resistance_marker: "NDM-1".to_string(),
            receptor_profile: probe_profile,
        };
        let probe_b = BacterialStrainSpec {
            name: "probe_b".to_string(),
            species: "Klebsiella pneumoniae".to_string(),
            resistance_marker: "KPC".to_string(),
            receptor_profile: probe_profile,
        };
        let pa = pairwise_score(&probe_a, &phages[0], temperature_k, coeffs);
        let pb = pairwise_score(&probe_b, &phages[0], temperature_k, coeffs);
        (pa.lysis_potential_score - pb.lysis_potential_score).abs()
    };

    PhageMatchingPanel {
        rows,
        best_by_strain,
        resistance_independence_probe_abs_delta,
        mean_best_lysis_score,
    }
}

pub fn default_phage_host_matching_panel(temperature_k: f64) -> PhageMatchingPanel {
    evaluate_phage_host_matching_panel(
        &default_strain_specs(),
        &default_phage_specs(),
        temperature_k,
        PhageMatchingCoefficients::default(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_size_matches_cartesian_product() {
        let strains = default_strain_specs();
        let phages = default_phage_specs();
        let panel = evaluate_phage_host_matching_panel(
            &strains,
            &phages,
            310.15,
            PhageMatchingCoefficients::default(),
        );
        assert_eq!(panel.rows.len(), strains.len() * phages.len());
    }

    #[test]
    fn resistance_independence_probe_is_zero() {
        let panel = default_phage_host_matching_panel(310.15);
        assert!(panel.resistance_independence_probe_abs_delta <= 1.0e-12);
    }

    #[test]
    fn ndm_klebsiella_prefers_omp_targeting_phage() {
        let panel = default_phage_host_matching_panel(310.15);
        let best = panel
            .best_by_strain
            .iter()
            .find(|b| b.strain_name == "kp_ndm1_clinical")
            .expect("ndm strain");
        assert_eq!(best.best_phage_name, "phi_kp_omp");
    }

    #[test]
    fn json_ingest_parses_minimal_payload() {
        let data = r#"{
            "strains": [{
                "name":"s1",
                "species":"Klebsiella pneumoniae",
                "resistance_marker":"NDM-1",
                "receptors":{"lamb":0.1,"ompk35":0.7,"ompk36":0.8,"fhua":0.2,"lps_core":0.9,"type_iv_pilus":0.1}
            }],
            "phages": [{
                "name":"p1",
                "family":"myoviridae_like",
                "primary_receptor":"OmpK36",
                "secondary_receptor":"OmpK35",
                "secondary_weight":0.5,
                "ionic_contact_count":2.0,
                "hbond_contact_count":3.0,
                "hydrophobic_area_a2":300.0,
                "conformational_entropy_penalty":2.0,
                "host_takeover_efficiency":0.8
            }]
        }"#;
        let strains = load_strains_from_json_str(data).expect("strains");
        let phages = load_phages_from_json_str(data).expect("phages");
        assert_eq!(strains.len(), 1);
        assert_eq!(phages.len(), 1);
    }

    #[test]
    fn csv_ingest_parses_minimal_payload() {
        let strains_csv = "name,species,resistance_marker,lamb,ompk35,ompk36,fhua,lps_core,type_iv_pilus\n\
s1,Klebsiella pneumoniae,NDM-1,0.1,0.7,0.8,0.2,0.9,0.1\n";
        let phages_csv = "name,family,primary_receptor,secondary_receptor,secondary_weight,ionic_contact_count,hbond_contact_count,hydrophobic_area_a2,conformational_entropy_penalty,host_takeover_efficiency\n\
p1,myoviridae_like,OmpK36,OmpK35,0.5,2.0,3.0,300.0,2.0,0.8\n";
        let strains = load_strains_from_csv_str(strains_csv).expect("strains");
        let phages = load_phages_from_csv_str(phages_csv).expect("phages");
        assert_eq!(strains.len(), 1);
        assert_eq!(phages.len(), 1);
    }
}
