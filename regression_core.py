"""Import, analysis and export for Regression 2 (no agreement-analysis branch)."""
from __future__ import annotations

import io
import json
import math
import re
import warnings
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score
from threadpoolctl import threadpool_limits

import retained_methods as original

BUILD = "Regression 2 / 2.1"
BOOTSTRAP_SEED = 20260722


def read_upload(name, payload):
    source = io.BytesIO(payload)
    if name.lower().endswith(".csv"):
        frames = {"CSV": pd.read_csv(source, float_precision="round_trip")}
    elif name.lower().endswith(".xlsx"):
        frames = pd.read_excel(source, sheet_name=None, engine="openpyxl")
    else:
        raise ValueError("Use .xlsx or .csv. Save legacy .xls workbooks as .xlsx first.")
    result = {}
    for name, frame in frames.items():
        # Preserve source position before the source app's deterministic sort.
        frame.columns = [str(c).strip() for c in frame.columns]
        if any(c.startswith("__r2_") for c in frame):
            raise ValueError("Column names starting with __r2_ are reserved by Regression 2.")
        frame["__r2_source_row"] = np.arange(2, len(frame) + 2)
        result[name] = original.canonical_rows(frame)
    return result


def _key(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def find_column(columns, candidates):
    for candidate in candidates:
        matches = [c for c in columns if _key(c) == _key(candidate)]
        if len(matches) == 1:
            return matches[0]
    return None


def detect_metadata(data):
    columns = [c for c in data if not c.startswith("__r2_")]
    detected = {
        "sample": find_column(columns, ["bloodSampleId", "Blood Sample ID", "sample_id", "Sample ID"]),
        "donor": find_column(columns, ["Donor", "donor_id", "donorid", "subject_id", "patient_id"]),
        "specimen": find_column(columns, ["specimen_type", "sample_type", "Type", "Specimen"]),
        "replicate": find_column(columns, ["replicate_number", "Replicate", "rep"]),
        "flag": find_column(columns, ["global_flag", "Global flag"]),
        "batch": find_column(columns, ["batch_id", "Batch ID"]),
    }
    return detected


def has_structured_sample_ids(data, mapping):
    """Only strip replicate suffixes when the entire nonblank ID layout supports it.

    A number after '-Hg-' can identify a sample or level, not a replicate.
    Mixed legacy/nonlegacy IDs must therefore use their full text consistently.
    """
    col = mapping.get("sample")
    if not col:
        return False
    values = _ids(data[col])
    present = values.notna()
    if not present.any():
        return False
    parsed = pd.DataFrame([original.parse_blood_sample_id(x) for x in values]).set_axis(data.index)
    if not parsed.loc[present, "parse_ok"].all():
        return False
    known_types = {value for group in original.DEFAULT_ANALYSIS_GROUPS.values() for value in group}
    known_layout = parsed.loc[present, "specimen_type"].isin(known_types).all()
    if mapping.get("replicate"):
        explicit = pd.to_numeric(data[mapping["replicate"]], errors="coerce")
        both = present & explicit.notna() & parsed["replicate_number"].notna()
        if both.any() and not explicit.loc[both].eq(parsed.loc[both, "replicate_number"]).all():
            return False
        return bool(known_layout or (both.sum() == present.sum()))
    return bool(known_layout)


def detect_pairs(columns):
    columns = [str(c) for c in columns if not str(c).startswith("__r2_")]
    rows = []
    for measured in columns:
        if re.search(r"(?:^|[_\s])(ref|reference|sysmex)(?:$|[_\s])", measured, re.I):
            continue
        base = re.sub(r"_\d+$", "", measured)
        candidates = [f"{measured}_ref", f"{base}_ref", f"{measured}_reference",
                      f"{base}_reference", f"ref_{measured}", f"reference_{measured}",
                      f"{base}_sysmex", f"sysmex_{base}"]
        ref = find_column(columns, candidates)
        if ref and ref != measured:
            cfg = original.ANALYTES.get(base.upper(), {})
            normal = cfg.get("normal", (np.nan, np.nan))
            rows.append({"Run": True, "Measured (Y)": measured, "Reference (X)": ref,
                         "Unit": cfg.get("unit", ""), "Normal low": normal[0], "Normal high": normal[1]})
    return pd.DataFrame(rows, columns=["Run", "Measured (Y)", "Reference (X)", "Unit", "Normal low", "Normal high"])


def validate_pairs(table, columns):
    selected = table.loc[table["Run"].fillna(False).eq(True)].copy()
    if selected.empty:
        raise ValueError("Enable at least one measured/reference pair.")
    used = set()
    for _, row in selected.iterrows():
        y, x = row["Measured (Y)"], row["Reference (X)"]
        if y not in columns or x not in columns:
            raise ValueError("Choose existing columns for both X and Y in every enabled row.")
        if y == x:
            raise ValueError(f"{y}: X and Y must be distinct, independently measured columns.")
        if y in used:
            raise ValueError(f"{y}: select each measured column only once.")
        used.add(y)
    return selected.reset_index(drop=True)


def _ids(series):
    return series.map(original.identifier).replace({"": None, "nan": None, "None": None, "<NA>": None})


def prepare_metadata(data, mapping, identity="auto"):
    if data.empty:
        raise ValueError("The selected worksheet has no data rows.")
    out = data.copy().reset_index(drop=True)
    if "__r2_source_row" not in out:
        out["__r2_source_row"] = np.arange(2, len(out) + 2)
    for col in mapping.values():
        if col is not None and col not in data:
            raise ValueError(f"Metadata column not found: {col}")
    # Raw batch identifiers are preserved, but removals always target unique rows.
    out["source_row"] = out["__r2_source_row"]
    out["source_batch_id"] = out[mapping["batch"]] if mapping.get("batch") else ""
    out["batch_id"] = [f"row_{i:08d}" for i in out["source_row"]]
    out["bloodSampleId"] = out[mapping["sample"]] if mapping.get("sample") else ""
    parsed = pd.DataFrame([original.parse_blood_sample_id(x) for x in out["bloodSampleId"]])
    structured_ids = has_structured_sample_ids(data, mapping)
    explicit = _ids(data[mapping["donor"]]).reset_index(drop=True) if mapping.get("donor") else None
    notes = []
    if identity == "auto":
        identity = "explicit" if explicit is not None else ("prefix" if structured_ids else "full_sample" if mapping.get("sample") else "row")
    if identity == "explicit":
        if explicit is None:
            raise ValueError("Choose a donor/sample column for explicit identity.")
        if structured_ids:
            _, errors = original.evaluate_donor_identity_mapping(parsed, explicit, "explicit_column")
            if errors:
                raise ValueError(" ".join(errors))
        out["donor"] = explicit
    elif identity == "full_sample":
        if not mapping.get("sample"):
            raise ValueError("Choose a sample-ID column to use full sample identities.")
        out["donor"] = _ids(data[mapping["sample"]]).reset_index(drop=True)
        notes.append("Repeated readings are grouped by the full sample ID. Numeric suffixes are preserved. Sample IDs are treated as independent units for confidence intervals; select an explicit donor column if multiple IDs belong to the same donor.")
    elif identity in {"prefix", "token"}:
        if not mapping.get("sample"):
            raise ValueError("Choose a parsable sample-ID column for this identity rule.")
        out["donor"] = parsed["prefix_donor_key" if identity == "prefix" else "parsed_donor_token"]
        issues, errors = original.evaluate_donor_identity_mapping(parsed, explicit, "prefix_plus_token" if identity == "prefix" else "parsed_token")
        if errors or (issues and identity == "token"):
            raise ValueError(" ".join(issues + errors) + " Choose an explicit donor column or collection prefix + D-token.")
        notes.extend(issues)
    elif identity == "row":
        out["donor"] = [f"Sample_{i}" for i in out["source_row"]]
        notes.append("Each input row is one independent unit. Select a donor/sample column when rows are repeated measurements.")
    else:
        raise ValueError("Unknown donor identity rule.")
    if mapping.get("specimen"):
        out["specimen_type"] = data[mapping["specimen"]].map(original.normalize_specimen_label).reset_index(drop=True)
    elif mapping.get("sample") and (identity in {"prefix", "token"} or (identity == "explicit" and structured_ids)):
        out["specimen_type"] = parsed["specimen_type"]
    else:
        out["specimen_type"] = "All"
    fallback = out.groupby(["donor", "specimen_type"], dropna=False).cumcount() + 1
    use_parsed_replicates = identity in {"prefix", "token"} or (identity == "explicit" and structured_ids)
    rep = pd.to_numeric(data[mapping["replicate"]], errors="coerce").reset_index(drop=True) if mapping.get("replicate") else parsed["replicate_number"] if use_parsed_replicates else pd.Series(np.nan, index=out.index)
    out["replicate_number"] = rep.where(rep.notna(), fallback)
    out["parse_ok"] = out[["donor", "specimen_type"]].notna().all(axis=1)
    if mapping.get("flag"):
        out["global_flag"] = data[mapping["flag"]].reset_index(drop=True)
    else:
        out["global_flag"] = False
        notes.append("No global-flag column selected: all rows are treated as FALSE.")
    out["global_flag_bool"] = out["global_flag"].map(original.normalize_bool)
    out["donor_identity_rule"] = identity
    out["parsed_donor_token"] = parsed["parsed_donor_token"]
    out["prefix_donor_key"] = parsed["prefix_donor_key"]
    return out, notes


def default_groups(types):
    rows, assigned = [], set()
    for name, labels in original.DEFAULT_ANALYSIS_GROUPS.items():
        found = [t for t in labels if t in types]
        if found:
            rows.append({"Run": True, "Group": name, "Specimen types": ", ".join(found)})
            assigned.update(found)
    rows += [{"Run": True, "Group": t, "Specimen types": t} for t in types if t not in assigned]
    return pd.DataFrame(rows, columns=["Run", "Group", "Specimen types"])


def validate_groups(table, types):
    groups = {}
    for _, row in table.iterrows():
        if not row.get("Run", False) or pd.isna(row.get("Run")):
            continue
        name = str(row["Group"]).strip()
        labels = [x.strip() for x in str(row["Specimen types"]).split(",") if x.strip()]
        if not name or name in groups or not labels or set(labels) - set(types):
            raise ValueError("Each enabled group needs a unique name and existing specimen types separated by commas.")
        groups[name] = labels
    if not groups:
        raise ValueError("Enable at least one regression group.")
    return groups


def normality(values, group, analyte, stage):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    row = {"analysis_group": group, "analyte": analyte, "stage": stage, "n": len(values),
           "shapiro_W": np.nan, "shapiro_p": np.nan, "alpha": 0.05}
    if len(values) < 3:
        row["interpretation"] = "Unavailable: fewer than 3 values"
    elif np.ptp(values) == 0:
        row["interpretation"] = "Unavailable: constant values"
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w, p = stats.shapiro(values)
        row.update(shapiro_W=float(w), shapiro_p=float(p), interpretation="Normality rejected" if p < 0.05 else "Normality not rejected")
        if len(values) > 5000:
            row["interpretation"] += "; p-value approximate for n > 5000"
    return row


def regression_metrics(x, y, group, analyte, reference, unit, n_boot):
    n = len(x)
    row = {"analysis_group": group, "analyte": analyte, "reference_column": reference,
           "unit": unit, "n_donors": n}
    names = ["huber_slope", "huber_slope_ci95_low", "huber_slope_ci95_high", "huber_intercept",
             "huber_intercept_ci95_low", "huber_intercept_ci95_high", "pearson_r", "pearson_r_ci95_low",
             "pearson_r_ci95_high", "pearson_p", "pearson_r_squared", "huber_fit_r_squared",
             "huber_fit_rmse", "r_squared_vs_identity", "rmse_vs_identity"]
    row.update({key: np.nan for key in names})
    row.update(bootstrap_requested=n_boot, bootstrap_successful=0, fit_warning_count=0,
               slope_intercept_ci_method="Donor-pair percentile bootstrap (2.5%, 97.5%)",
               correlation_ci_method="Fisher-z", bootstrap_seed=BOOTSTRAP_SEED, status="OK")
    if n < 2 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        row["status"] = "Unavailable: fewer than 2 donor pairs or constant X/Y"
        return row
    r, p = stats.pearsonr(x, y)
    lo, hi = original.pearson_ci(r, n)
    row.update(pearson_r=float(r), pearson_p=float(p), pearson_r_ci95_low=lo,
               pearson_r_ci95_high=hi, pearson_r_squared=float(r * r),
               r_squared_vs_identity=float(r2_score(x, y)), rmse_vs_identity=float(np.sqrt(mean_squared_error(x, y))))
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            slope, intercept = original.huber_regression(x, y)
            slo, shi, ilo, ihi, successful = original.bootstrap_huber_ci(x, y, n_boot=n_boot)
        fitted = slope * x + intercept
        row.update(huber_slope=slope, huber_intercept=intercept,
                   huber_equation=f"Y = {slope:.8g} X {intercept:+.8g}",
                   huber_slope_ci95_low=slo, huber_slope_ci95_high=shi,
                   huber_intercept_ci95_low=ilo, huber_intercept_ci95_high=ihi,
                   bootstrap_successful=successful, fit_warning_count=len(caught),
                   huber_fit_r_squared=float(r2_score(y, fitted)),
                   huber_fit_rmse=float(np.sqrt(mean_squared_error(y, fitted))))
        notes = []
        if n < 4:
            notes.append("95% CIs unavailable: at least 4 independent donor pairs required")
        elif successful < n_boot:
            notes.append(f"{n_boot-successful} bootstrap resamples failed or had constant X")
        if successful < 0.9 * n_boot and n >= 4:
            notes.append("CI based on fewer than 90% of requested resamples; review")
        if caught:
            notes.append(f"{len(caught)} fitting warnings; inspect convergence before interpreting intervals")
        if notes:
            row["status"] = "; ".join(notes)
    except Exception as exc:
        row["status"] = f"Huber fit unavailable: {exc}"
    return row


@threadpool_limits.wrap(limits=1)
def analyze(data, prepared, pairs, groups, donors=None, n_boot=5000, progress=None):
    """Run every pair independently so shared reference columns cannot collide."""
    if not 1 <= n_boot <= 20000:
        raise ValueError("Bootstrap resamples must be between 1 and 20,000.")
    frames = {name: [] for name in ["Regression Results", "Donor Means", "Outliers", "ESD Diagnostics", "Normality", "Missing Pairs"]}
    keep = prepared["parse_ok"] & prepared["specimen_type"].isin({v for types in groups.values() for v in types})
    if donors is not None:
        keep &= prepared["donor"].isin(donors)
    selection_excluded = prepared.loc[~keep].copy()
    selection_excluded["exclusion_reason"] = np.where(~selection_excluded["parse_ok"], "Missing/unparsed donor or specimen", "Donor/specimen not selected")
    selected = prepared.loc[keep].copy()
    if selected.empty:
        raise ValueError("No rows remain after donor/specimen selection.")
    flag_excluded = selected.loc[~selected["global_flag_bool"].eq(False)].copy()
    flag_excluded["exclusion_reason"] = np.where(flag_excluded["global_flag_bool"].eq(True), "global_flag = TRUE", "global_flag missing or unrecognized")
    eligible = selected.loc[selected["global_flag_bool"].eq(False)].copy()
    if eligible.empty:
        raise ValueError("No rows have global_flag = FALSE after selection.")
    plot_data = []
    total = len(groups) * len(pairs)
    counter = 0
    for group, types in groups.items():
        group_rows = eligible.loc[eligible["specimen_type"].isin(types)].copy()
        for _, pair in pairs.iterrows():
            ycol, xcol = pair["Measured (Y)"], pair["Reference (X)"]
            unit = "" if pd.isna(pair.get("Unit", "")) else str(pair.get("Unit", ""))
            work = group_rows.copy()
            # Internal pair names isolate all source columns from metadata aliases.
            work["__r2_x"] = original.finite_numeric(data.loc[work.index, xcol])
            work["__r2_y"] = original.finite_numeric(data.loc[work.index, ycol])
            cfg = {ycol: {"mhs": "__r2_y", "ref": "__r2_x"}}
            missing = work.loc[work[["__r2_x", "__r2_y"]].isna().any(axis=1)].copy()
            if not missing.empty:
                missing["analysis_group"], missing["analyte"] = group, ycol
                missing["exclusion_reason"] = "Missing/non-numeric/non-finite X or Y; paired row excluded for this analyte"
                frames["Missing Pairs"].append(missing)
            complete = work.dropna(subset=["__r2_x", "__r2_y"])
            frames["Normality"].append(pd.DataFrame([normality(complete["__r2_y"]-complete["__r2_x"], group, ycol, "Raw paired differences before ESD")]))
            if len(complete) >= 20 and (complete[["__r2_x", "__r2_y"]] <= -10001).any().any():
                row = regression_metrics(np.array([]), np.array([]), group, ycol, xcol, unit, n_boot)
                row["status"] = "Not analyzed: inherited ESD log1p(value + 10000) requires values > -10001"
                frames["Regression Results"].append(pd.DataFrame([row]))
            else:
                removed, diagnostic = original.detect_esd_outliers(work, group, analytes=cfg)
                if not removed.empty:
                    removed = removed.merge(work[["batch_id", "source_row", "source_batch_id"]], on="batch_id", how="left", validate="many_to_one")
                    removed["mhs_column"], removed["reference_column"] = ycol, xcol
                frames["Outliers"].append(removed)
                diagnostic["mhs_column"], diagnostic["reference_column"] = ycol, xcol
                frames["ESD Diagnostics"].append(diagnostic)
                if work.empty:
                    means = pd.DataFrame(columns=["Donor", "__r2_x", "__r2_y", "Replicates", f"{ycol}_n"])
                else:
                    means, counts = original.build_donor_means(work, removed, analytes=cfg)
                    means = means.merge(counts, on="Donor", how="left").dropna(subset=["__r2_x", "__r2_y"])
                x, y = means["__r2_x"].to_numpy(float), means["__r2_y"].to_numpy(float)
                row = regression_metrics(x, y, group, ycol, xcol, unit, n_boot)
                row.update(n_complete_replicates=len(complete), n_esd_removed=len(removed), n_missing_pairs=len(missing))
                frames["Regression Results"].append(pd.DataFrame([row]))
                means = means.rename(columns={"__r2_x": "reference_mean", "__r2_y": "measured_mean", f"{ycol}_n": "retained_replicates"})
                means["analysis_group"], means["analyte"], means["reference_column"] = group, ycol, xcol
                frames["Donor Means"].append(means)
                residuals = y - (row["huber_slope"] * x + row["huber_intercept"])
                frames["Normality"].append(pd.DataFrame([normality(residuals, group, ycol, "Donor-mean Huber fit residuals after ESD")]))
                plot_data.append({"row": row, "x": x, "y": y, "residuals": residuals,
                                  "normal": [pair.get("Normal low", np.nan), pair.get("Normal high", np.nan)]})
            counter += 1
            if progress:
                progress(counter / total, f"{group}: {ycol} ({counter}/{total})")
    tables = {}
    for name, items in frames.items():
        nonempty = [frame for frame in items if not frame.empty]
        tables[name] = pd.concat(nonempty, ignore_index=True) if nonempty else (items[0].iloc[:0].copy() if items else pd.DataFrame())
    tables.update({"global flag exclusions": flag_excluded, "Selection Exclusions": selection_excluded,
                   "Input Mapping": pairs, "Donor Identity": prepared[["source_row", "source_batch_id", "bloodSampleId", "donor", "specimen_type", "replicate_number", "parse_ok", "donor_identity_rule", "parsed_donor_token", "prefix_donor_key"]]})
    analysis_unit = "sample ID" if prepared["donor_identity_rule"].eq("full_sample").all() else "input row" if prepared["donor_identity_rule"].eq("row").all() else "donor/sample identity"
    tables["Regression Results"]["analysis_unit"] = analysis_unit
    tables["Regression Results"]["n_independent_units"] = tables["Regression Results"]["n_donors"]
    tables["Methods"] = pd.DataFrame([
        ("App", BUILD), ("Regression direction", "X = reference; Y = measured; fit on retained paired donor means"),
        ("Grouping unit", f"{analysis_unit}; identity rule: {', '.join(prepared['donor_identity_rule'].unique())}. n_independent_units is the number of paired means. n_donors and Donor Means are retained legacy field/sheet names and do not verify donor independence."),
        ("Global flag", "Only FALSE retained; TRUE and missing/unrecognized values excluded. If column absent, all rows treated as FALSE."),
        ("Outliers", "Original generalized ESD per analyte on raw linked replicate rows, before donor averaging; original log1p(value+10000), SD/CV/mixed branches retained; minimum 20 rows per tested subset."),
        ("ESD details", "Alpha 0.05; maximum fraction 0.05 (0.08 in the original two-bin special branch). No removal if the tested subset has fewer than 20 rows."),
        ("Huber estimator", "Original sklearn HuberRegressor() defaults: epsilon=1.35, alpha=0.0001, max_iter=100, tol=1e-5."),
        ("Slope/intercept 95% CI", f"Original paired donor bootstrap with {n_boot} resamples, seed {BOOTSTRAP_SEED}; percentile 2.5/97.5. ESD and donor averaging are fixed before resampling."),
        ("Pearson R 95% CI", "Original Fisher-z interval using n independent donor means and SE=1/sqrt(n-3); at least 4 donors. Near-perfect R clipped to ±0.999999999 as in source."),
        ("Normality", "Shapiro-Wilk alpha 0.05 on raw paired differences and final Huber residuals; diagnostic only. The original regression ESD branch does not switch to Grubbs/MAD."),
        ("R-squared definitions", "pearson_r_squared = R²; huber_fit_r_squared = 1-SSE(Y, fitted Y)/SST(Y). Legacy r_squared_vs_identity = sklearn.r2_score(X,Y). These are distinct."),
        ("RMSE definitions", "huber_fit_rmse compares Y to fitted Y; legacy rmse_vs_identity compares Y to X."),
        ("Independence", "Each donor contributes one paired mean within a group. Row-as-donor mode assumes independent input rows. Overlapping groups are separate analyses, not independent cohorts."),
        ("CI interpretation", "Bootstrap intervals are conditional on selected rows and ESD removals; uncertainty from outlier selection is not re-estimated. Fitting warnings and successful resample counts are reported."),
        ("Missing results", "Unavailable estimates stay blank/NaN with a status. Correlation and robust regression do not alone establish interchangeability."),
        ("Source", "User-supplied Interchangeability-App-main (8).zip, active app.py embedded regression engine 2026-09-27-proxima-v7-validation."),
        ("Documentation", "https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.HuberRegressor.html ; https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats._result_classes.PearsonRResult.confidence_interval.html"),
    ], columns=["Setting", "Description"])
    return {"tables": tables, "plots": plot_data}


def _number(value):
    return f"{value:.5g}" if np.isfinite(value) else "unavailable"


def draw_regression(item, ax, show_ranges=True):
    row, x, y = item["row"], item["x"], item["y"]
    ax.scatter(x, y, s=27, color="#257c8b", alpha=.85, edgecolor="white", linewidth=.35)
    if len(x):
        left, right = float(np.min(x)), float(np.max(x))
        pad = max((right-left)*.05, abs(left)*.005, .01)
        grid = np.linspace(left-pad, right+pad, 200)
        ax.plot(grid, grid, "--", color="#9ca3af", lw=1, label="Identity")
        if np.isfinite(row["huber_slope"]):
            ax.plot(grid, row["huber_slope"]*grid+row["huber_intercept"], color="#203e63", lw=1.8, label="Huber fit")
        if show_ranges:
            for val in item["normal"]:
                if pd.notna(val) and left-pad <= float(val) <= right+pad:
                    ax.axvline(float(val), color="#bd6c6c", lw=.9, alpha=.65)
        ax.legend(fontsize=8, loc="best", frameon=False)
    ax.set_title(f"{row['analyte']} · {row['analysis_group']}", fontsize=11, pad=9)
    unit = f" ({row['unit']})" if row["unit"] else ""
    ax.set_xlabel(f"{row['reference_column']}{unit}", fontsize=9)
    ax.set_ylabel(f"{row['analyte']}{unit}", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=.13)
    lines = [f"n = {row['n_donors']}   |   R = {_number(row['pearson_r'])}  [95% CI {_number(row['pearson_r_ci95_low'])}, {_number(row['pearson_r_ci95_high'])}]",
             f"Slope = {_number(row['huber_slope'])}  [95% CI {_number(row['huber_slope_ci95_low'])}, {_number(row['huber_slope_ci95_high'])}]",
             f"Intercept = {_number(row['huber_intercept'])}  [95% CI {_number(row['huber_intercept_ci95_low'])}, {_number(row['huber_intercept_ci95_high'])}]"]
    ax.text(0, -.25, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=8.5, linespacing=1.6)


def regression_figure(item, show_ranges=True):
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    fig.subplots_adjust(left=.13, right=.97, top=.9, bottom=.32)
    draw_regression(item, ax, show_ranges)
    return fig


def diagnostic_figure(item):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), layout="constrained")
    res = item["residuals"]
    if len(res) and np.isfinite(res).all():
        fitted = item["y"] - res
        axes[0].scatter(fitted, res, color="#257c8b", s=24)
        axes[0].axhline(0, ls="--", lw=1, color="#9ca3af")
        if len(res) >= 3 and np.ptp(res) > 0:
            stats.probplot(res, dist="norm", plot=axes[1])
        else:
            axes[1].text(.5, .5, "Q-Q plot unavailable", ha="center", transform=axes[1].transAxes)
    axes[0].set(xlabel="Fitted Y", ylabel="Y − fitted Y", title="Huber fit residuals")
    axes[1].set_title("Normal Q-Q plot")
    fig.suptitle(f"{item['row']['analyte']} · {item['row']['analysis_group']}", fontsize=11)
    return fig


def plot_bytes(fig, fmt="png"):
    stream = io.BytesIO()
    fig.savefig(stream, format=fmt, dpi=300, facecolor="white", bbox_inches="tight")
    return stream.getvalue()


def export_plots(result, show_ranges=True):
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for i, item in enumerate(result["plots"]):
            stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{i+1:03d}_{item['row']['analysis_group']}_{item['row']['analyte']}")
            for kind, fig in [("regression", regression_figure(item, show_ranges)), ("diagnostics", diagnostic_figure(item))]:
                for fmt in ["png", "pdf"]:
                    archive.writestr(f"{stem}_{kind}.{fmt}", plot_bytes(fig, fmt))
                plt.close(fig)
        for page in range(0, len(result["plots"]), 6):
            items = result["plots"][page:page+6]
            if not items:
                continue
            nrows = math.ceil(len(items)/2)
            fig, axes = plt.subplots(nrows, 2, figsize=(14, 6.3*nrows), squeeze=False)
            fig.subplots_adjust(left=.08, right=.97, top=.94, bottom=.30/nrows+.025, wspace=.29, hspace=.72)
            for ax, item in zip(axes.flat, items):
                draw_regression(item, ax, show_ranges)
            for ax in list(axes.flat)[len(items):]:
                ax.set_visible(False)
            for fmt in ["png", "pdf"]:
                archive.writestr(f"Regression_2_panel_{page//6+1}.{fmt}", plot_bytes(fig, fmt))
            plt.close(fig)
    return output.getvalue()


def export_workbook(result):
    """Preserve the original app's Python/XlsxWriter export stack."""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine="xlsxwriter", engine_kwargs={"options": {"strings_to_formulas": False, "strings_to_urls": False}}) as writer:
        book = writer.book
        header = book.add_format({"bold": True, "bg_color": "#203E63", "font_color": "white", "text_wrap": True, "valign": "top"})
        number = book.add_format({"num_format": "0.00000"})
        pvalue = book.add_format({"num_format": "0.000E+00"})
        wrap = book.add_format({"text_wrap": True, "valign": "top"})
        for name, frame in result["tables"].items():
            frame = frame.copy()
            # Internal columns are not part of the user's source schema.
            frame = frame[[c for c in frame if not str(c).startswith("__r2_")]]
            if len(frame.columns) == 0:
                frame = pd.DataFrame({"Status": ["No rows"]})
            frame.to_excel(writer, sheet_name=name, index=False)
            sheet = writer.sheets[name]
            sheet.freeze_panes(1, min(2, len(frame.columns)))
            sheet.autofilter(0, 0, max(1, len(frame)), len(frame.columns)-1)
            sheet.set_row(0, 44)
            for j, col in enumerate(frame):
                sheet.write(0, j, col, header)
                width = min(36, max(16, len(str(col))*.9))
                fmt = pvalue if col in {"pearson_p", "shapiro_p"} else number if pd.api.types.is_float_dtype(frame[col]) else None
                sheet.set_column(j, j, width, fmt)
            if name in {"Methods", "Input Audit"}:
                sheet.set_column(0, 0, 27, wrap)
                sheet.set_column(1, 1, 115, wrap)
                for i in range(1, len(frame)+1):
                    sheet.set_row(i, 46)
    return output.getvalue()
