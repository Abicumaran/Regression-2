from __future__ import annotations

import hashlib
import io
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from regression_core import (BUILD, analyze, default_groups, detect_metadata, detect_pairs,
                             diagnostic_figure, export_plots, export_workbook, prepare_metadata,
                             read_upload, regression_figure, validate_groups, validate_pairs)

st.set_page_config(page_title="Regression 2", page_icon="📈", layout="wide")
st.title("Regression 2")
st.caption("Huber regression · 95% confidence intervals for slope, intercept and Pearson’s R")

with st.expander("Analysis method"):
    st.markdown("""
    The analysis follows the regression branch of your original app: global-flag filtering,
    generalized ESD on raw paired replicates, analyte-specific removal, donor averaging,
    then Huber regression with **X = reference** and **Y = measured**.

    Slope and intercept use the original donor-pair percentile bootstrap. Pearson’s R uses
    the original Fisher-z interval. All three 95% confidence intervals are always requested.
    Shapiro–Wilk results and residual/Q-Q plots are diagnostics; they do not change the
    original regression estimator or ESD rule. At least four independent donor pairs are
    required for the confidence intervals. Missing estimates remain blank with an explanation.
    """)

sample_path = Path(__file__).with_name("sample_input_template.csv")
if sample_path.exists():
    st.download_button("Download example input (synthetic CSV)", sample_path.read_bytes(),
                       "Regression_2_example.csv", "text/csv")
uploaded = st.file_uploader("1. Upload an Excel workbook or CSV", type=["xlsx", "csv"])
demo = st.checkbox("Try the synthetic example", value=False) if uploaded is None else False
if uploaded is None and not demo:
    st.info("Upload data with measured columns and their matched reference columns, such as PLT_3 and PLT_ref. Other names can be mapped below.")
    st.stop()

name = uploaded.name if uploaded else sample_path.name
payload = uploaded.getvalue() if uploaded else sample_path.read_bytes()
digest = hashlib.sha256(payload).hexdigest()
try:
    sheets = read_upload(name, payload)
except Exception as exc:
    st.error(str(exc))
    st.stop()
sheet = st.selectbox("Worksheet", list(sheets), key=f"sheet_{digest}")
data = sheets[sheet]
cols = [c for c in data if not c.startswith("__r2_")]
st.caption(f"{len(data):,} rows · {len(cols):,} columns")
with st.expander("Input preview"):
    st.dataframe(data[cols].head(50), hide_index=True, use_container_width=True)

st.subheader("2. Data detection")
detected = detect_metadata(data)
mapping = {}
with st.expander("Review metadata mappings", expanded=False):
    options = ["<none>"] + cols
    labels = {"sample": "Sample ID", "donor": "Donor / independent sample ID",
              "specimen": "Specimen type", "replicate": "Replicate number",
              "flag": "Global flag", "batch": "Source batch / row ID"}
    boxes = st.columns(3)
    for i, (field, label) in enumerate(labels.items()):
        value = boxes[i % 3].selectbox(label, options, index=options.index(detected[field]) if detected[field] else 0,
                                      key=f"meta_{field}_{digest}_{sheet}")
        mapping[field] = None if value == "<none>" else value
    identities = {"Automatic (detect sample-ID format)": "auto",
                  "Full sample ID (keep text unchanged)": "full_sample",
                  "Explicit donor / sample column": "explicit", "Collection prefix + D-token": "prefix",
                  "Parsed D-token only": "token", "Each row is an independent sample": "row"}
    identity_label = st.selectbox("Donor identity rule", list(identities))
try:
    prepared, notices = prepare_metadata(data, mapping, identities[identity_label])
except Exception as exc:
    st.error(str(exc))
    st.stop()
for notice in notices:
    st.info(notice)
a, b, c = st.columns(3)
a.metric("Independent donors / samples", prepared.loc[prepared.parse_ok, "donor"].nunique())
b.metric("Specimen types", prepared.loc[prepared.parse_ok, "specimen_type"].nunique())
c.metric("Rows with global flag FALSE", int(prepared["global_flag_bool"].eq(False).sum()))
if (~prepared.parse_ok).any():
    st.warning(f"{int((~prepared.parse_ok).sum())} rows have missing/unparsed identities and will be excluded.")
with st.expander("Donor identity audit"):
    st.dataframe(prepared[["source_row", "source_batch_id", "bloodSampleId", "donor", "specimen_type", "replicate_number", "parse_ok"]], hide_index=True, use_container_width=True)
types = sorted(prepared.loc[prepared.parse_ok, "specimen_type"].astype(str).unique())
donors = sorted(prepared.loc[prepared.parse_ok, "donor"].astype(str).unique())
selected_donors = st.multiselect("Donors / samples to include", donors, default=donors, key=f"donors_{digest}_{sheet}")

st.subheader("3. Measured and reference columns")
st.caption("Common *_ref pairs and separate model columns are detected automatically. Add or edit rows for other datasets. X and Y must use compatible units.")
pairs = detect_pairs(cols)
if pairs.empty:
    pairs = pd.DataFrame([{"Run": True, "Measured (Y)": None, "Reference (X)": None,
                           "Unit": "", "Normal low": None, "Normal high": None}])
edited_pairs = st.data_editor(pairs, num_rows="dynamic", hide_index=True, use_container_width=True,
    key=f"pairs_{digest}_{sheet}", column_config={
        "Run": st.column_config.CheckboxColumn(default=True),
        "Measured (Y)": st.column_config.SelectboxColumn(options=cols, required=True),
        "Reference (X)": st.column_config.SelectboxColumn(options=cols, required=True),
        "Normal low": st.column_config.NumberColumn(help="Optional plot marker only; never filters data."),
        "Normal high": st.column_config.NumberColumn(help="Optional plot marker only; never filters data."),
    })
st.caption("Normal-range markers retain the source app’s defaults for recognized analytes. Edit or clear them for another population. They do not determine residual normality.")

st.subheader("4. Regression groups and confidence intervals")
group_table = st.data_editor(default_groups(types), num_rows="dynamic", hide_index=True, use_container_width=True,
    key=f"groups_{digest}_{sheet}", column_config={"Run": st.column_config.CheckboxColumn(default=True)})
st.caption("Mc and Mc2 are combined by default, as in the source app. Other specimen types retain separate groups.")
a, b = st.columns(2)
n_boot = a.number_input("Donor-level bootstrap resamples", min_value=200, max_value=20000, value=5000, step=200,
                         help="Original bootstrap method and fixed seed 20260722. More resamples take longer.")
show_ranges = b.checkbox("Show normal-range markers on plots", value=True)
st.caption("Automatic generalized ESD is enabled. For tested subsets with fewer than 20 complete replicate rows, the original method performs no ESD removal.")

signature = hashlib.sha256(json.dumps({"input": digest, "sheet": sheet, "mapping": mapping,
    "identity": identities[identity_label], "donors": selected_donors, "pairs": edited_pairs.to_json(),
    "groups": group_table.to_json(), "n_boot": n_boot, "ranges": show_ranges}, sort_keys=True).encode()).hexdigest()
run = st.button("Run Regression 2", type="primary", use_container_width=True)
if run:
    try:
        selected_pairs = validate_pairs(edited_pairs, cols)
        groups = validate_groups(group_table, types)
        if not selected_donors:
            raise ValueError("Select at least one donor/sample.")
        progress = st.progress(0, text="Starting regression analysis…")
        with st.spinner("Fitting Huber models and calculating confidence intervals…"):
            result = analyze(data, prepared, selected_pairs, groups, selected_donors, int(n_boot),
                             progress=lambda value, text: progress.progress(value, text=text))
            result["tables"]["Methods"].loc[len(result["tables"]["Methods"])] = ["Input", f"{name}; worksheet {sheet}; SHA256 {digest}"]
            workbook = export_workbook(result)
        st.session_state["regression_2_output"] = {"signature": signature, "result": result, "workbook": workbook}
        st.session_state.pop("regression_2_plots", None)
        progress.empty()
    except Exception as exc:
        st.error(f"Analysis could not finish: {exc}")
        st.stop()

saved = st.session_state.get("regression_2_output")
if not saved:
    st.stop()
if saved["signature"] != signature:
    st.info("Settings or input changed. Run Regression 2 again to refresh the results and downloads.")
    st.stop()
result = saved["result"]
st.divider()
st.subheader("Regression results")
metrics = result["tables"]["Regression Results"]
if not metrics.empty:
    statuses = metrics.loc[metrics.status.ne("OK"), ["analysis_group", "analyte", "status"]]
    if not statuses.empty:
        st.warning("Some results need review or have unavailable intervals. See status details below.")
        st.dataframe(statuses, hide_index=True, use_container_width=True)
    primary = ["analysis_group", "analyte", "analysis_unit", "n_independent_units", "huber_slope", "huber_slope_ci95_low", "huber_slope_ci95_high",
               "huber_intercept", "huber_intercept_ci95_low", "huber_intercept_ci95_high", "pearson_r", "pearson_r_ci95_low", "pearson_r_ci95_high"]
    st.dataframe(metrics[primary], hide_index=True, use_container_width=True)
st.download_button("Download Regression 2 Excel report", saved["workbook"], "Regression_2_Results.xlsx",
                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
st.download_button("Download all regression metrics (CSV)", metrics.to_csv(index=False).encode(),
                   "Regression_2_Metrics.csv", "text/csv")
tabs = st.tabs(["Plots", "Normality", "Outliers and exclusions", "All metrics", "Methods"])
with tabs[0]:
    if result["plots"]:
        choice = st.selectbox("Plot", range(len(result["plots"])),
            format_func=lambda i: f"{result['plots'][i]['row']['analysis_group']} — {result['plots'][i]['row']['analyte']}")
        item = result["plots"][choice]
        fig = regression_figure(item, show_ranges)
        st.pyplot(fig, use_container_width=False)
        plt.close(fig)
        with st.expander("Residual and Q-Q plots"):
            fig = diagnostic_figure(item)
            st.pyplot(fig)
            plt.close(fig)
        if st.button("Prepare all plots (PNG + PDF, individual + panels)"):
            with st.spinner("Preparing plots…"):
                st.session_state["regression_2_plots"] = export_plots(result, show_ranges)
        if "regression_2_plots" in st.session_state:
            st.download_button("Download plots ZIP", st.session_state["regression_2_plots"], "Regression_2_Plots.zip", "application/zip")
with tabs[1]:
    st.caption("Shapiro–Wilk is diagnostic. A non-significant p-value does not prove normality.")
    st.dataframe(result["tables"]["Normality"], use_container_width=True, hide_index=True)
with tabs[2]:
    for key in ["Outliers", "ESD Diagnostics", "global flag exclusions", "Selection Exclusions", "Missing Pairs", "Donor Means"]:
        with st.expander(key):
            st.dataframe(result["tables"][key], use_container_width=True, hide_index=True)
with tabs[3]:
    st.dataframe(metrics, use_container_width=True, hide_index=True)
    st.caption("R² of the Huber fit, R² of the correlation, and the source app’s R² against identity have separate labels.")
with tabs[4]:
    st.dataframe(result["tables"]["Methods"], use_container_width=True, hide_index=True)
st.caption(BUILD)
