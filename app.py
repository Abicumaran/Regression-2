from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from regression_core import (default_groups, detect_metadata, prepare_metadata,
                             read_upload, validate_groups, validate_pairs)
from analysis_job import AnalysisJob
from results_ui import show_results

from ui_controls import pair_controls, group_controls, display_table, table_styles

st.set_page_config(page_title="Regression 2", page_icon="📈", layout="wide")
table_styles()
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

    Bias is **mean(measured − reference)**. Bias (%) divides that signed bias by
    the mean reference value and multiplies by 100; absolute bias (%) is its absolute value.
    The main bias table uses retained paired readings after filtering and ESD. A separate
    donor-mean table uses the same paired means as the regression. These calculations
    compare measured values with references, not predictions from the fitted Huber line.
    """)

sample_path = Path(__file__).with_name("sample_input_template.csv")
if sample_path.exists():
    st.download_button("Download example input (synthetic CSV)", sample_path.read_bytes(),
                       "Regression_2_example.csv", "text/csv")
uploaded = st.file_uploader("1. Upload an Excel workbook or CSV", type=["xlsx", "csv"])
demo = st.checkbox("Try the synthetic example", value=False) if uploaded is None else False
if uploaded is None and not demo:
    st.info("Upload data with measured columns and their matched reference columns, such as PLT_3 and PLT_ref. Other names can be mapped below.")
    show_results()
    st.stop()

name = uploaded.name if uploaded else sample_path.name
payload = uploaded.getvalue() if uploaded else sample_path.read_bytes()
digest = hashlib.sha256(payload).hexdigest()
try:
    sheets = read_upload(name, payload)
except Exception as exc:
    st.error(str(exc))
    show_results()
    st.stop()
sheet = st.selectbox("Worksheet", list(sheets), key=f"sheet_{digest}")
data = sheets[sheet]
cols = [c for c in data if not c.startswith("__r2_")]
st.caption(f"{len(data):,} rows · {len(cols):,} columns")
with st.expander("Input preview"):
    display_table(data[cols].head(50), hide_index=True, use_container_width=True)

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
    show_results()
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
    display_table(prepared[["source_row", "source_batch_id", "bloodSampleId", "donor", "specimen_type", "replicate_number", "parse_ok"]], hide_index=True, use_container_width=True)
types = sorted(prepared.loc[prepared.parse_ok, "specimen_type"].astype(str).unique())
donors = sorted(prepared.loc[prepared.parse_ok, "donor"].astype(str).unique())
selected_donors = st.multiselect("Donors / samples to include", donors, default=donors, key=f"donors_{digest}_{sheet}")

st.subheader("3. Measured and reference columns")
edited_pairs = pair_controls(cols, f"pairs_{digest}_{sheet}")
st.caption("Normal-range markers retain the source app’s defaults for recognized analytes. Edit or clear them for another population. They do not determine residual normality.")

st.subheader("4. Regression groups and confidence intervals")
group_table = group_controls(default_groups(types), types, f"groups_{digest}_{sheet}")
st.caption("Mc and Mc2 are combined by default, as in the source app. Other specimen types retain separate groups.")
a, b = st.columns(2)
n_boot = a.number_input("Donor-level bootstrap resamples", min_value=200, max_value=20000, value=5000, step=200,
                         help="Original bootstrap method and fixed seed 20260722. More resamples take longer.")
show_ranges = b.checkbox("Show normal-range markers on plots", value=True)
st.caption("Automatic generalized ESD is enabled. For tested subsets with fewer than 20 complete replicate rows, the original method performs no ESD removal.")

signature = hashlib.sha256(json.dumps({"input": digest, "sheet": sheet, "mapping": mapping,
    "identity": identities[identity_label], "donors": selected_donors, "pairs": edited_pairs.to_json(),
    "groups": group_table.to_json(), "n_boot": n_boot, "ranges": show_ranges}, sort_keys=True).encode()).hexdigest()
job = st.session_state.get("regression_2_job")
busy = job is not None and not job.snapshot()["done"]
run = st.button("Run Regression 2", type="primary", use_container_width=True, disabled=busy)
if run and not busy:
    try:
        selected_pairs = validate_pairs(edited_pairs, cols)
        groups = validate_groups(group_table, types)
        if not selected_donors:
            raise ValueError("Select at least one donor/sample.")
        job = AnalysisJob(data, prepared, selected_pairs, groups, selected_donors, int(n_boot),
                          signature=signature, source=f"{name}; worksheet {sheet}; SHA256 {digest}",
                          show_ranges=show_ranges)
        st.session_state["regression_2_job"] = job
        job.start()
        st.rerun()
    except Exception as exc:
        st.error(f"Analysis could not start: {exc}")

show_results(signature)
