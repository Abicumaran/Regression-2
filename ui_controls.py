"""Plain Streamlit controls and HTML tables, without the lazy-loaded data grid."""
import pandas as pd
import streamlit as st

from regression_core import detect_pairs, find_column

ANALYTE_ORDER = (
    "RBC_3", "WBC_2", "PLT_3", "HCT", "HGB", "MCV_3", "RDW_3",
    "MCH", "MCHC", "NEUT_2", "LYPH_2", "MXD_2",
)
PAIR_COLUMNS = ["Run", "Measured (Y)", "Reference (X)", "Unit", "Normal low", "Normal high"]


def analyte_choices(columns):
    """Keep the requested panel order, then append other detected measurements."""
    detected = detect_pairs(columns).to_dict("records")
    choices, used = [], set()
    for label in ANALYTE_ORDER:
        measured = find_column(columns, [label])
        pair = next((row.copy() for row in detected if row["Measured (Y)"] == measured), None)
        if pair is None:
            pair = dict(zip(PAIR_COLUMNS, [False, measured, None, "", None, None]))
        if measured is not None:
            used.add(measured)
        choices.append((label, pair))
    choices.extend((row["Measured (Y)"], row) for row in detected if row["Measured (Y)"] not in used)
    return choices


def _optional_number(value):
    return None if value is None or pd.isna(value) else float(value)


def pair_controls(columns, key):
    choices = analyte_choices(columns)
    st.caption("Uncheck an analyte to exclude it. Missing columns stay unselected; expand its mapping to choose your column names. Additional detected analytes follow the requested panel.")
    extra = st.number_input("Additional custom analytes", min_value=0, max_value=50, value=0, key=f"{key}_extras")
    choices += [(f"Custom analyte {i + 1}", dict(zip(PAIR_COLUMNS, [False, None, None, "", None, None]))) for i in range(int(extra))]
    options = [None] + list(columns)
    rows = []
    for i, (label, pair) in enumerate(choices):
        prefix = f"{key}_{i}"
        toggle, settings = st.columns([1, 4])
        enabled = toggle.checkbox(f"Include {label}", value=bool(pair["Run"]), key=f"{prefix}_run")
        current_y = st.session_state.get(f"{prefix}_y", pair["Measured (Y)"])
        current_x = st.session_state.get(f"{prefix}_x", pair["Reference (X)"])
        summary = f"Y: {current_y} · X: {current_x}" if current_y is not None and current_x is not None else "Choose measured and reference columns"
        with settings.expander(f"{label} — {summary}"):
            left, right = st.columns(2)
            measured = left.selectbox("Measured (Y)", options, index=options.index(pair["Measured (Y)"]),
                                      format_func=lambda x: "Choose a column" if x is None else x, key=f"{prefix}_y")
            reference = right.selectbox("Reference (X)", options, index=options.index(pair["Reference (X)"]),
                                        format_func=lambda x: "Choose a column" if x is None else x, key=f"{prefix}_x")
            unit = st.text_input("Unit", value=pair["Unit"], key=f"{prefix}_unit")
            left, right = st.columns(2)
            low = left.number_input("Normal low (optional)", value=_optional_number(pair["Normal low"]), key=f"{prefix}_low")
            high = right.number_input("Normal high (optional)", value=_optional_number(pair["Normal high"]), key=f"{prefix}_high")
        rows.append(dict(zip(PAIR_COLUMNS, [enabled, measured, reference, unit, low, high])))
    return pd.DataFrame(rows, columns=PAIR_COLUMNS)


def group_controls(defaults, types, key):
    rows = defaults.to_dict("records")
    extra = st.number_input("Additional custom groups", min_value=0, max_value=20, value=0, key=f"{key}_extras")
    rows += [{"Run": False, "Group": f"Custom group {i + 1}", "Specimen types": ""} for i in range(int(extra))]
    edited = []
    for i, row in enumerate(rows):
        prefix = f"{key}_{i}"
        enabled = st.checkbox(f"Include group: {row['Group']}", value=row["Run"], key=f"{prefix}_run")
        with st.expander(f"{row['Group']} — group settings"):
            name = st.text_input("Group name", value=row["Group"], key=f"{prefix}_name")
            selected = st.multiselect("Specimen types", types,
                default=[x.strip() for x in row["Specimen types"].split(",") if x.strip()], key=f"{prefix}_types")
        edited.append({"Run": enabled, "Group": name, "Specimen types": ", ".join(selected)})
    return pd.DataFrame(edited, columns=["Run", "Group", "Specimen types"])


def display_table(frame, **_):
    """Escape uploaded content and bound the on-screen preview; exports stay complete."""
    if frame.empty:
        st.caption("No rows to display.")
        return
    html = frame.head(200).to_html(index=False, escape=True, na_rep="", border=0,
                                   float_format=lambda value: f"{value:.6g}", classes="r2-table")
    st.markdown('<div class="r2-table-wrap" tabindex="0" role="region" aria-label="Data table">' + html + '</div>', unsafe_allow_html=True)
    if len(frame) > 200:
        st.caption(f"Showing the first 200 of {len(frame):,} rows. Download the report for all rows.")


def table_styles():
    st.markdown("""<style>
    .r2-table-wrap {overflow:auto; max-height:480px; margin:0.5rem 0 1rem;
      border:1px solid #8885; border-radius:8px;}
    .r2-table {border-collapse:collapse; width:100%; font-size:0.875rem;}
    .r2-table th,.r2-table td {padding:0.6rem 0.8rem; text-align:left;
      border-bottom:1px solid #8883; white-space:nowrap;}
    .r2-table th {font-weight:600;}
    </style>""", unsafe_allow_html=True)
