"""Persistent results panel, separate from input validation and analysis execution."""
import matplotlib.pyplot as plt
import streamlit as st

from regression_core import BUILD, diagnostic_figure, export_plots, export_workbook, regression_figure
from ui_controls import display_table


def render_results(saved, current_signature):
    show_ranges = saved.get("show_ranges", True)
    result = saved["result"]
    st.divider()
    st.subheader("Regression results")
    st.success("Analysis complete — results saved for this session.")
    st.caption(f"Completed: {saved.get('completed', 'earlier in this session')} · {saved.get('source', '')}")
    st.caption("Download the Excel report or CSV to keep a permanent copy on your computer. A browser reload, closed tab, or server restart may end this session.")
    if current_signature != saved["signature"]:
        st.info("Settings or input changed. These are the last completed results, using the previous settings. Run Regression 2 again to replace them.")
    if saved.get("export_error"):
        st.warning(f"Analysis is saved, but the Excel report could not be prepared: {saved['export_error']}. You can still download the metrics CSV.")
    if saved.get("export_error") and st.button("Retry Excel report"):
        try:
            saved["workbook"] = export_workbook(result)
            saved["export_error"] = None
            st.rerun()
        except Exception as exc:
            st.error(f"Excel export failed: {exc}. The calculated results remain available.")
    if saved.get("workbook") is not None:
        st.download_button("Download Regression 2 Excel report", saved["workbook"], "Regression_2_Results.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary", on_click="ignore")
    st.download_button("Download all regression metrics (CSV)", result["tables"]["Regression Results"].to_csv(index=False).encode(),
                       "Regression_2_Metrics.csv", "text/csv", on_click="ignore")
    metrics = result["tables"]["Regression Results"]
    if not metrics.empty:
        if "bias" in metrics:
            st.subheader("Bias — retained paired readings")
            st.caption("Positive bias means measured values are higher than the reference. Bias (%) = 100 × mean(measured − reference) / mean(reference). Absolute bias (%) = |bias (%)|. Percentages are unavailable when mean reference is zero.")
            bias_columns = ["analysis_group", "analyte", "unit", "n_bias_pairs", "mean_true", "bias", "bias_percent", "abs_bias_percent", "bias_status"]
            bias_labels = {"analysis_group": "Group", "analyte": "Analyte", "unit": "Unit",
                           "n_bias_pairs": "Paired readings", "mean_true": "Mean reference", "bias": "Bias",
                           "bias_percent": "Bias (%)", "abs_bias_percent": "Absolute bias (%)", "bias_status": "Status"}
            display_table(metrics[bias_columns].rename(columns=bias_labels))
            with st.expander("Bias on donor means — same population as the regression"):
                st.caption("Each retained donor/sample contributes one paired mean. This can differ from reading-level bias when replicate counts differ.")
                donor_columns = ["analysis_group", "analyte", "unit"] + [f"donor_mean_{key}" for key in bias_columns[3:]]
                donor_labels = {**bias_labels, **{f"donor_mean_{key}": value for key, value in bias_labels.items()}}
                donor_labels["donor_mean_n_bias_pairs"] = "Paired donor means"
                display_table(metrics[donor_columns].rename(columns=donor_labels))
        else:
            st.info("Run Regression 2 again to calculate bias for these saved results.")
        st.subheader("Regression estimates")
        statuses = metrics.loc[metrics.status.ne("OK"), ["analysis_group", "analyte", "status"]]
        if not statuses.empty:
            st.warning("Some results need review or have unavailable intervals. See status details below.")
            display_table(statuses, hide_index=True, use_container_width=True)
        primary = ["analysis_group", "analyte", "analysis_unit", "n_independent_units", "huber_slope", "huber_slope_ci95_low", "huber_slope_ci95_high",
                   "huber_intercept", "huber_intercept_ci95_low", "huber_intercept_ci95_high", "pearson_r", "pearson_r_ci95_low", "pearson_r_ci95_high"]
        display_table(metrics[primary], hide_index=True, use_container_width=True)
    else:
        st.warning("Analysis finished with no regression rows. Review your analyte, donor, and group selections.")
    tabs = st.tabs(["Plots", "Normality", "Outliers and exclusions", "All metrics", "Methods"])
    with tabs[0]:
        if result["plots"]:
            choice = st.selectbox("Plot", range(len(result["plots"])),
                format_func=lambda i: f"{result['plots'][i]['row']['analysis_group']} — {result['plots'][i]['row']['analyte']}")
            item = result["plots"][choice]
            try:
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
                    st.download_button("Download plots ZIP", st.session_state["regression_2_plots"], "Regression_2_Plots.zip", "application/zip", on_click="ignore")
            except Exception as exc:
                st.error(f"Plot could not be prepared: {exc}. Your results and report downloads remain available.")
    with tabs[1]:
        st.caption("Shapiro–Wilk is diagnostic. A non-significant p-value does not prove normality.")
        display_table(result["tables"]["Normality"], use_container_width=True, hide_index=True)
    with tabs[2]:
        for key in ["Outliers", "ESD Diagnostics", "global flag exclusions", "Selection Exclusions", "Missing Pairs", "Donor Means"]:
            with st.expander(key):
                display_table(result["tables"][key], use_container_width=True, hide_index=True)
    with tabs[3]:
        display_table(metrics, use_container_width=True, hide_index=True)
        st.caption("R² of the Huber fit, R² of the correlation, and the source app’s R² against identity have separate labels.")
    with tabs[4]:
        display_table(result["tables"]["Methods"], use_container_width=True, hide_index=True)
    st.caption(BUILD)


def show_results(current_signature=None):
    job = st.session_state.get("regression_2_job")
    polling = job is not None and not job.handled

    @st.fragment(run_every=1 if polling else None)
    def panel():
        if job is not None:
            state = job.snapshot()
            if state["output"] is not None:
                if st.session_state.get("regression_2_output", {}).get("result") is not state["output"]["result"]:
                    st.session_state.pop("regression_2_plots", None)
                st.session_state["regression_2_output"] = state["output"]
            if not state["done"]:
                st.progress(state["progress"], text=state["message"])
                st.info("Analysis is running. Results will appear here automatically. You can adjust controls without interrupting this run; it uses the settings submitted when you clicked Run.")
            elif not job.handled:
                job.handled = True
                st.rerun()  # Stop polling and re-enable the Run button.
            if state["error"]:
                st.error(f"Analysis could not finish: {state['error']}. Any previous completed results are kept below.")
        saved = st.session_state.get("regression_2_output")
        if saved:
            render_results(saved, current_signature)

    panel()
