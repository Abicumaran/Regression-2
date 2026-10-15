"""Regression 2: retained statistical methods from the supplied app's active V9 engine.

Functions below are copied verbatim for numerical parity. The separate
Bland-Altman engine is intentionally not included.
"""
from __future__ import annotations
import json
import re
import warnings
from typing import Mapping, Sequence
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t
from sklearn.linear_model import HuberRegressor

ANALYTES: dict[str, dict[str, object]] = {
    "RBC":   {"mhs": "RBC",     "ref": "RBC_ref",   "normal": (3.85, 5.65),   "unit": "10^12/L"},
    "WBC":   {"mhs": "WBC_2",   "ref": "WBC_ref",   "normal": (3.60, 10.50),  "unit": "10^9/L"},
    "NEUT":  {"mhs": "NEUT_2",  "ref": "NEUT_ref",  "normal": (1.50, 7.70),   "unit": "10^9/L"},
    "LYMPH": {"mhs": "LYMPH_2", "ref": "LYMPH_ref", "normal": (1.10, 4.00),   "unit": "10^9/L"},
    "MXD":   {"mhs": "MXD_2",   "ref": "MXD_ref",   "normal": (0.10, 1.60),   "unit": "10^9/L"},
    "PLT":   {"mhs": "PLT_3",     "ref": "PLT_ref",   "normal": (160.0, 370.0), "unit": "10^9/L"},
    "MCV":   {"mhs": "MCV_3",     "ref": "MCV_ref",   "normal": (80.0, 101.0),  "unit": "fL"},
    "RDW":   {"mhs": "RDW_3",     "ref": "RDW_ref",   "normal": (11.0, 16.0),   "unit": "%"},
    "MCH":   {"mhs": "MCH",     "ref": "MCH_ref",   "normal": (27.0, 34.0),   "unit": "pg"},
    "MCHC":  {"mhs": "MCHC",    "ref": "MCHC_ref",  "normal": (320.0, 360.0), "unit": "g/L"},
    "HGB":   {"mhs": "HGB",     "ref": "HGB_ref",   "normal": (118.0, 172.0), "unit": "g/L"},
    "HCT":   {"mhs": "HCT",     "ref": "HCT_ref",   "normal": (35.0, 49.0),   "unit": "percentage points"},
}

DEFAULT_ANALYSIS_GROUPS: dict[str, list[str]] = {
    "Capillary_Mc_Mc2": ["Mc", "Mc2"],
    "Venous_MK2": ["MK2"],
    "Venous_MK3": ["MK3"],
    "Direct_Capillary_Dir": ["Dir"],
    "TASSO_Tas": ["Tas"],
}


def identifier(value):
    if pd.isna(value):
        return None
    if isinstance(value, (int, np.integer)):
        return str(value)
    if isinstance(value, (float, np.floating)) and np.isfinite(value) and value.is_integer():
        return str(int(value))
    return str(value).strip()

def canonical_rows(df):
    """Unique positional index and deterministic rows; never silently deduplicate."""
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    if out.columns.duplicated().any():
        raise ValueError('Duplicate column names after trimming spaces.')
    keys = [c for c in ['Condition', 'Interferent', 'Level', 'bloodSampleId',
                        'Blood Sample ID', 'sample_id', 'donor', 'specimen_type',
                        'Day', 'Device', 'deviceId', 'Replicate', 'batch_id'] if c in out]
    for c in keys:
        if c not in ['Replicate']:
            out[c] = out[c].map(identifier)
    if keys:
        # A string sort keeps IDs stable across mixed Excel and CSV dtypes.
        order = out[keys].fillna('').astype(str).sort_values(keys, kind='stable').index
        # Duplicate source indices are common after concatenating workbooks.
        tmp = out.reset_index(drop=True)
        order = tmp[keys].fillna('').astype(str).sort_values(keys, kind='stable').index
        out = tmp.loc[order]
    return out.reset_index(drop=True)

def finite_numeric(values):
    return pd.to_numeric(values, errors='coerce').replace([np.inf, -np.inf], np.nan)

def normalize_bool(value) -> bool | None:
    if pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None

def normalize_specimen_label(value: object) -> object:
    """Normalize common specimen aliases while preserving unknown labels."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    text = str(value).strip().strip("-_")
    if not text:
        return np.nan
    aliases = {
        "mc": "Mc", "mc2": "Mc2", "mk2": "MK2", "mk3": "MK3",
        "dir": "Dir", "direct": "Dir", "tas": "Tas", "tasso": "Tas",
    }
    return aliases.get(text.lower(), text)

def parse_blood_sample_id(value: object) -> dict[str, object]:
    """Parse donor/specimen/replicate metadata from a blood-sample ID.

    The parsed donor token starts with ``D`` and ends immediately before the
    first hyphen. A collection/date prefix is also retained so the app can
    detect cases such as ``140725D03`` versus ``130726D03``. The final donor
    identity is selected by the app: explicit Donor column, parsed D-token, or
    prefix + D-token.
    """
    text = "" if value is None else str(value).strip()
    donor_match = re.search(r"(?P<donor>D\d+[A-Za-z0-9]*)-", text, flags=re.IGNORECASE)
    if not donor_match:
        return {
            "donor": np.nan, "parsed_donor_token": np.nan,
            "sample_prefix": np.nan, "prefix_donor_key": np.nan,
            "specimen_type": np.nan, "replicate_number": np.nan,
            "parse_ok": False,
        }

    donor = donor_match.group("donor")
    # Standardize the leading D while preserving suffix case exactly as entered.
    donor = "D" + donor[1:]
    prefix = text[:donor_match.start()].strip()
    prefix_key = f"{prefix}{donor}" if prefix else donor
    tail = text[donor_match.end():].strip()
    parts = [part.strip() for part in tail.split("-") if part.strip()]
    if not parts:
        return {
            "donor": donor, "parsed_donor_token": donor,
            "sample_prefix": prefix, "prefix_donor_key": prefix_key,
            "specimen_type": np.nan, "replicate_number": np.nan,
            "parse_ok": False,
        }

    replicate = np.nan
    if len(parts) >= 2 and re.fullmatch(r"\d+", parts[-1]):
        replicate = int(parts[-1])
        specimen = "-".join(parts[:-1]).strip()
    else:
        specimen = "-".join(parts).strip()
    specimen = normalize_specimen_label(specimen)

    return {
        "donor": donor,
        "parsed_donor_token": donor,
        "sample_prefix": prefix,
        "prefix_donor_key": prefix_key,
        "specimen_type": specimen,
        "replicate_number": replicate,
        "parse_ok": bool(pd.notna(specimen)),
    }

def generalized_esd_with_batch(x, batch_ids, alpha: float = 0.05, max_outliers_frac: float = 0.05) -> dict:
    x = np.asarray(x, dtype=float)
    batch_ids = np.asarray(batch_ids)
    mask = ~np.isnan(x)
    x = x[mask]
    batch_ids = batch_ids[mask]
    n_total = len(x)
    if n_total < 20:
        raise ValueError("ESD requires at least 20 samples")
    h = int(np.floor(max_outliers_frac * n_total))
    if h < 1:
        return {"n_outliers": 0, "outlier_batch_ids": [], "test_statistics": [], "critical_values": []}

    remaining = x.copy()
    remaining_ids = batch_ids.copy()
    esd_stats, lambdas, removed_ids = [], [], []
    for i in range(1, h + 1):
        mean = remaining.mean()
        sd = remaining.std(ddof=1)
        if sd == 0:
            break
        deviations = np.abs(remaining - mean) / sd
        max_idx = int(np.argmax(deviations))
        esd_i = float(deviations[max_idx])
        nu = n_total - i - 1
        p = 1 - alpha / (2 * (n_total - i + 1))
        t_crit = t.ppf(p, nu)
        lambda_i = ((n_total - i) * t_crit) / np.sqrt((nu + t_crit**2) * (n_total - i + 1))
        esd_stats.append(esd_i)
        lambdas.append(float(lambda_i))
        removed_ids.append(remaining_ids[max_idx])
        remaining = np.delete(remaining, max_idx)
        remaining_ids = np.delete(remaining_ids, max_idx)

    n_outliers = 0
    for i, (esd_i, lambda_i) in enumerate(zip(esd_stats, lambdas), start=1):
        if esd_i > lambda_i:
            n_outliers = i
    return {
        "n_outliers": n_outliers,
        "outlier_batch_ids": [str(v) for v in removed_ids[:n_outliers]],
        "test_statistics": esd_stats,
        "critical_values": lambdas,
    }

def sd_cv_pattern(y_true, y_pred, n_bins: int = 5, min_per_bin: int = 5, eps: float = 1e-8):
    y_true = np.asarray(y_true, float).reshape(-1)
    y_pred = np.asarray(y_pred, float).reshape(-1)
    if y_true.size != y_pred.size:
        raise ValueError("y_true and y_pred must have the same length")
    if not (np.isfinite(y_true).all() and np.isfinite(y_pred).all()):
        raise ValueError("There is some NaN/inf")
    d = y_true - y_pred
    m = 0.5 * (y_true + y_pred)
    edges = np.unique(np.quantile(m, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 4:
        raise ValueError("Too few unique values to create bins")
    bin_id = np.digitize(m, edges[1:-1], right=True)
    rows = []
    for k in range(edges.size - 1):
        mask = bin_id == k
        if mask.sum() < min_per_bin:
            continue
        mk = m[mask]
        dk = d[mask]
        sd = dk.std(ddof=1)
        mean_m = mk.mean()
        cv = sd / max(abs(mean_m), eps)
        rows.append({"bin": k, "range": (edges[k], edges[k + 1]), "n": int(mask.sum()),
                     "mean_level_m": float(mean_m), "sd_diff": float(sd), "cv_local": float(cv)})
    tbl = pd.DataFrame(rows)
    min_valid_bins = max(2, n_bins // 2)
    if tbl.shape[0] < min_valid_bins:
        raise ValueError(f"Very few valid bins (<{min_valid_bins})")
    sd_ratio = tbl["sd_diff"].max() / max(tbl["sd_diff"].min(), eps)
    cv_ratio = tbl["cv_local"].max() / max(tbl["cv_local"].min(), eps)
    corr_sd = np.corrcoef(tbl["mean_level_m"], tbl["sd_diff"])[0, 1]
    corr_cv = np.corrcoef(tbl["mean_level_m"], tbl["cv_local"])[0, 1]
    if sd_ratio <= 2.0 and abs(corr_sd) < 0.3 and corr_cv < -0.3:
        label, lbl = "SD constant (additive)", "SD constant"
    elif cv_ratio <= 2.0 and abs(corr_cv) < 0.3 and corr_sd > 0.3:
        label, lbl = "CV constant (percent)", "CV constant"
    else:
        label, lbl = "Mixed or uncertain (review by low/high sections)", "Mixed"
    return tbl, {
        "sd_ratio_maxmin": float(sd_ratio), "cv_ratio_maxmin": float(cv_ratio),
        "corr(level, sd)": float(corr_sd), "corr(level, cv)": float(corr_cv),
        "diagnosis": label, "lbl": lbl, "edges_used": edges.tolist(),
    }, lbl

def change_point_sd_by_level(y_true, y_pred, n_bins: int = 5, min_per_bin: int = 5,
                             robust: bool = True, eps: float = 1e-8):
    y_true = np.asarray(y_true, float).reshape(-1)
    y_pred = np.asarray(y_pred, float).reshape(-1)
    d = y_true - y_pred
    m = 0.5 * (y_true + y_pred)
    edges = np.unique(np.quantile(m, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 4:
        raise ValueError("Not enough unique m values to create bins")
    bin_id = np.digitize(m, edges[1:-1], right=True)
    rows = []
    for k in range(edges.size - 1):
        mask = bin_id == k
        nk = int(mask.sum())
        if nk < min_per_bin:
            continue
        dk = d[mask]
        mk = m[mask]
        if robust:
            med = np.median(dk)
            mad = np.median(np.abs(dk - med))
            scale = 1.4826 * mad
        else:
            scale = dk.std(ddof=1)
        rows.append({"bin": k, "left": float(edges[k]), "right": float(edges[k + 1]),
                     "n": nk, "mean_level_m": float(mk.mean()), "scale": float(scale)})
    tbl = pd.DataFrame(rows).reset_index(drop=True)
    min_valid_bins = max(2, n_bins // 2)
    if tbl.shape[0] < min_valid_bins:
        raise ValueError(f"Too few valid bins (<{min_valid_bins})")
    log_scales = np.log(tbl["scale"].values + eps)
    jumps = np.abs(np.diff(log_scales))
    cut_pos = int(np.argmax(jumps))
    cutoff_m = float(tbl.loc[cut_pos, "right"])
    return cutoff_m, tbl, {
        "method": "robust MAD scale" if robust else "SD scale",
        "cut_between_bins": (int(tbl.loc[cut_pos, "bin"]), int(tbl.loc[cut_pos + 1, "bin"])),
        "jump_score_abs_log": float(jumps[cut_pos]), "cutoff_m": cutoff_m,
        "edges_used": edges.tolist(),
    }

def run_esd(diff, batch_ids, alpha: float = 0.05, max_outliers_frac: float = 0.05):
    if len(diff) < 20:
        return [], 0, {"n_outliers": 0, "outlier_batch_ids": [], "test_statistics": [], "critical_values": []}
    result = generalized_esd_with_batch(diff, batch_ids, alpha=alpha, max_outliers_frac=max_outliers_frac)
    return result["outlier_batch_ids"], result["n_outliers"], result

def run_binned_esd(df: pd.DataFrame, true_col: str, pred_col: str, batch_col: str = "batch_id",
                   n_bins: int = 10, max_outliers_frac: float = 0.05):
    work = df.copy()
    work["bin"] = pd.qcut(work[true_col], q=n_bins, duplicates="drop")
    outlier_ids, details = [], []
    for b in work["bin"].dropna().unique():
        sub = work[work["bin"] == b]
        if len(sub) < 20:
            continue
        diff = sub[pred_col] - sub[true_col]
        out_bin, _, result = run_esd(diff, sub[batch_col], max_outliers_frac=max_outliers_frac)
        outlier_ids.extend(out_bin)
        details.append({"bin": str(b), "result": result})
    return outlier_ids, len(outlier_ids), details

def compute_diff(df: pd.DataFrame, true_col: str, pred_col: str, method: str) -> pd.Series:
    if method == "SD":
        return df[pred_col] - df[true_col]
    if method == "CV":
        return (df[pred_col] - df[true_col]) / df[true_col] * 100
    raise ValueError("Unknown diff method")

def _validated_esd_single_analyte(work: pd.DataFrame, true_col: str, pred_col: str,
                                  n_bins: int = 2, max_outliers_bins_frac: float = 0.08) -> tuple[list[str], dict]:
    original_n = len(work)
    if original_n < 20:
        return [], {"error_model": "Not run: <20 replicate rows", "n_replicates": original_n,
                    "n_outliers": 0, "cutoff": np.nan, "esd_detail": "[]"}

    transformed = work[["batch_id", true_col, pred_col]].dropna().copy()
    if len(transformed) < 20:
        return [], {"error_model": "Not run: <20 complete replicate rows", "n_replicates": len(transformed),
                    "n_outliers": 0, "cutoff": np.nan, "esd_detail": "[]"}

    # Preserve the original validated transformation exactly.
    transformed[true_col] = np.log1p(transformed[true_col] + 10000)
    transformed[pred_col] = np.log1p(transformed[pred_col] + 10000)

    n_size = len(transformed)
    if n_size >= 300:
        n_size_bins, min_per_bin = 10, 30
    elif n_size >= 100:
        n_size_bins, min_per_bin = 5, 30
    elif n_size >= 50:
        n_size_bins, min_per_bin = 4, 15
    else:
        n_size_bins, min_per_bin = 4, 5

    try:
        _, summary, error_model = sd_cv_pattern(
            transformed[true_col], transformed[pred_col], n_bins=n_size_bins, min_per_bin=min_per_bin
        )
    except Exception as exc:
        summary, error_model = {"diagnosis": f"Fallback due to: {exc}"}, "Mixed"

    try:
        cutoff, _, cp_info = change_point_sd_by_level(
            transformed[true_col], transformed[pred_col], n_bins=n_size_bins,
            min_per_bin=min_per_bin, robust=True
        )
    except Exception as exc:
        cutoff, cp_info = float(transformed[true_col].median()), {"error": str(exc)}

    esd_detail: object = []
    if error_model == "SD constant":
        diff = compute_diff(transformed, true_col, pred_col, "SD")
        outlier_ids, _, esd_detail = run_esd(diff, transformed["batch_id"])
    elif error_model == "CV constant":
        diff = compute_diff(transformed, true_col, pred_col, "CV")
        outlier_ids, _, esd_detail = run_esd(diff, transformed["batch_id"])
    else:
        # This follows the original branching. With its +10000 log transform,
        # the near-zero branch is normally not entered.
        if (transformed[true_col] < 1).any():
            outlier_ids, _, esd_detail = run_binned_esd(
                transformed, true_col, pred_col, n_bins=n_bins,
                max_outliers_frac=max_outliers_bins_frac
            )
        else:
            low = transformed[transformed[true_col] <= cutoff]
            high = transformed[transformed[true_col] > cutoff]
            out_low, _, detail_low = run_esd(compute_diff(low, true_col, pred_col, "SD"), low["batch_id"])
            out_high, _, detail_high = run_esd(compute_diff(high, true_col, pred_col, "CV"), high["batch_id"])
            outlier_ids = list(out_low) + list(out_high)
            esd_detail = {"low": detail_low, "high": detail_high}

    diagnostics = {
        "error_model": error_model,
        "error_model_diagnosis": summary.get("diagnosis", ""),
        "n_replicates": len(transformed),
        "n_outliers": len(outlier_ids),
        "cutoff": cutoff,
        "change_point_info": json.dumps(cp_info, default=str),
        "esd_detail": json.dumps(esd_detail, default=str),
    }
    return [str(v) for v in outlier_ids], diagnostics

def detect_esd_outliers(group_df: pd.DataFrame, group_name: str,
                        analytes: Mapping[str, Mapping[str, object]] = ANALYTES,
                        n_bins: int = 2, max_outliers_bins_frac: float = 0.08) -> tuple[pd.DataFrame, pd.DataFrame]:
    outlier_rows, diagnostics_rows = [], []
    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        complete = group_df.dropna(subset=[mhs_col, ref_col]).copy()
        ids, diag = _validated_esd_single_analyte(
            complete, ref_col, mhs_col, n_bins=n_bins,
            max_outliers_bins_frac=max_outliers_bins_frac
        )
        diagnostics_rows.append({"analysis_group": group_name, "analyte": analyte,
                                 "mhs_column": mhs_col, "reference_column": ref_col, **diag})
        if ids:
            selected = group_df[group_df["batch_id"].astype(str).isin(ids)].copy()
            for _, row in selected.iterrows():
                ref = row.get(ref_col, np.nan)
                mhs = row.get(mhs_col, np.nan)
                residual = mhs - ref if pd.notna(mhs) and pd.notna(ref) else np.nan
                residual_pct = 100 * residual / ref if pd.notna(residual) and pd.notna(ref) and ref != 0 else np.nan
                outlier_rows.append({
                    "analysis_group": group_name,
                    "analyte": analyte,
                    "batch_id": str(row.get("batch_id", "")),
                    "bloodSampleId": row.get("bloodSampleId", ""),
                    "donor": row.get("donor", ""),
                    "specimen_type": row.get("specimen_type", ""),
                    "replicate_number": row.get("replicate_number", np.nan),
                    "mhs_column": mhs_col,
                    "reference_column": ref_col,
                    "mhs_value": mhs,
                    "reference_value": ref,
                    "residual_native": residual,
                    "residual_percent": residual_pct,
                    "exclusion_reason": "Validated generalized ESD; linked MHS/reference replicate removed for this analyte only",
                })
    outlier_columns = [
        "analysis_group", "analyte", "batch_id", "bloodSampleId", "donor",
        "specimen_type", "replicate_number", "mhs_column", "reference_column",
        "mhs_value", "reference_value", "residual_native", "residual_percent",
        "exclusion_reason",
    ]
    return pd.DataFrame(outlier_rows, columns=outlier_columns), pd.DataFrame(diagnostics_rows)

def build_donor_means(group_df: pd.DataFrame, outliers_df: pd.DataFrame,
                      analytes: Mapping[str, Mapping[str, object]] = ANALYTES) -> tuple[pd.DataFrame, pd.DataFrame]:
    donors = sorted(group_df["donor"].dropna().astype(str).unique())
    result = pd.DataFrame({"Donor": donors})
    total_counts = group_df.groupby("donor").size().rename("Replicates").reset_index().rename(columns={"donor": "Donor"})
    result = result.merge(total_counts, on="Donor", how="left")
    counts = pd.DataFrame({"Donor": donors})

    for analyte, cfg in analytes.items():
        mhs_col, ref_col = str(cfg["mhs"]), str(cfg["ref"])
        out_ids: set[str] = set()
        if not outliers_df.empty:
            out_ids = set(outliers_df.loc[outliers_df["analyte"] == analyte, "batch_id"].astype(str))
        sub = group_df[~group_df["batch_id"].astype(str).isin(out_ids)].dropna(subset=[mhs_col, ref_col]).copy()
        grouped = sub.groupby("donor")[[ref_col, mhs_col]].mean().reset_index().rename(columns={"donor": "Donor"})
        n = sub.groupby("donor").size().rename(f"{analyte}_n").reset_index().rename(columns={"donor": "Donor"})
        result = result.merge(grouped, on="Donor", how="left")
        counts = counts.merge(n, on="Donor", how="left")

    result.insert(0, "#", np.arange(1, len(result) + 1))
    ordered = ["#", "Donor"]
    for analyte, cfg in analytes.items():
        ordered.extend([str(cfg["ref"]), str(cfg["mhs"])])
    ordered.append("Replicates")
    result = result[ordered]
    return result, counts

def huber_regression(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    model = HuberRegressor()
    model.fit(np.asarray(x).reshape(-1, 1), np.asarray(y))
    return float(model.coef_[0]), float(model.intercept_)

def pearson_ci(r: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(r):
        return np.nan, np.nan
    r = float(np.clip(r, -0.999999999, 0.999999999))
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - z_crit * se)), float(np.tanh(z + z_crit * se))

def bootstrap_huber_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 5000,
                       alpha: float = 0.05, random_state: int = 20260722) -> tuple[float, float, float, float, int]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 4:
        return np.nan, np.nan, np.nan, np.nan, 0
    rng = np.random.default_rng(random_state)
    slopes, intercepts = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        xb, yb = x[idx], y[idx]
        if len(np.unique(xb)) < 2:
            continue
        try:
            slope, intercept = huber_regression(xb, yb)
            slopes.append(slope)
            intercepts.append(intercept)
        except Exception:
            continue
    if not slopes:
        return np.nan, np.nan, np.nan, np.nan, 0
    q = [100 * alpha / 2, 100 * (1 - alpha / 2)]
    slo, shi = np.percentile(slopes, q)
    ilo, ihi = np.percentile(intercepts, q)
    return float(slo), float(shi), float(ilo), float(ihi), len(slopes)

def evaluate_donor_identity_mapping(
    preliminary: pd.DataFrame,
    verification_donors: pd.Series | None,
    donor_identity_mode: str,
) -> tuple[list[str], list[str]]:
    """Return (issues needing human confirmation, invalid mapping errors).

    Multiple collection prefixes for the *same parsed D-token* do not by
    themselves make an explicit donor column ambiguous: the explicit labels
    can distinguish D03 from D03b.  Explicit labels also preserve separate preparations sharing a legacy prefix.
    Missing explicit labels remain a blocking error.
    """
    valid = preliminary.loc[preliminary["parse_ok"].fillna(False)].copy()
    if valid.empty:
        return [], ["No valid donor/sample IDs could be parsed."]

    mode = str(donor_identity_mode).strip().lower()
    explicit = None
    if verification_donors is not None:
        explicit = pd.Series(verification_donors).reset_index(drop=True).astype("string").str.strip()
        explicit = explicit.mask(explicit.isna() | explicit.str.lower().isin(["", "nan", "none", "null"]))
        valid["explicit"] = explicit.loc[valid.index].values

    warnings: list[str] = []
    errors: list[str] = []
    prefix_counts = valid.groupby("parsed_donor_token")["prefix_donor_key"].nunique()
    repeated_token = int((prefix_counts > 1).sum())
    if mode == "explicit_column":
        if explicit is None:
            errors.append("Choose the explicit Donor column or a different donor identity rule.")
        else:
            missing_count = int(valid["explicit"].isna().sum())
            if missing_count:
                errors.append(
                    f"The selected explicit Donor column is missing for {missing_count} parsed rows. "
                    "Fill the donor labels or choose a different identity rule; "
                    "do not silently fall back to a parsed D-token."
                )
            # Explicit specimen/donor labels are authoritative. They may intentionally
            # distinguish preparations sharing a legacy sample-ID prefix (e.g. D04b/D04m).
            # Never collapse those labels back to a parsed token.
    elif mode == "prefix_plus_token":
        if explicit is not None:
            splits = valid.dropna(subset=["explicit"]).groupby("explicit")["prefix_donor_key"].nunique()
            if int((splits > 1).sum()):
                warnings.append(
                    "One or more explicit donor labels span different collection prefixes; "
                    "prefix + D-token would split these labels into separate donors."
                )
        elif repeated_token:
            warnings.append(
                f"{repeated_token} D-token(s) recur under different collection prefixes. "
                "Confirm that those collection-specific identities should remain separate."
            )
    elif mode == "parsed_token":
        if repeated_token:
            warnings.append(
                f"{repeated_token} D-token(s) occur with multiple collection prefixes. "
                "Using the token alone will merge them into a single donor."
            )
        if explicit is not None:
            merged = valid.dropna(subset=["explicit"]).groupby("parsed_donor_token")["explicit"].nunique()
            if int((merged > 1).sum()):
                warnings.append(
                    "One or more parsed D-tokens map to different explicit Donor labels "
                    "(for example D03 versus D03b); token-only mode would merge them."
                )
    else:
        errors.append(f"Unsupported donor identity rule: {donor_identity_mode}")
    return warnings, errors
