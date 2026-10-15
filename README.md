# Regression 2

Live app: https://regression-v2.streamlit.app/

Version 2.1 replaces the editable data grids with checkboxes and ordinary column
selectors. Results and audits use escaped, scrollable HTML tables, avoiding the
grid JavaScript module that failed to load in some browsers. On-screen tables
preview at most 200 rows; downloadable reports contain the complete results.

The analyte checklist always starts in this order: **RBC_3, WBC_2, PLT_3, HCT,
HGB, MCV_3, RDW_3, MCH, MCHC, NEUT_2, LYPH_2, MXD_2**. Detected pairs are selected
automatically; missing pairs remain unchecked until the user maps real columns.
Uncheck any analyte to exclude it. Other detected analytes appear afterward, and
custom mappings and regression groups can still be added. Analysis, result rows,
and plot choices preserve the selected analyte order within each group.

After an update, reload an already-open tab to load the new app. A browser holding
the previous deployment may require Ctrl+Shift+R. The calculation methods are unchanged.

A regression-only Streamlit app adapted from the active regression engine in the
supplied `Interchangeability-App-main (8).zip`.

Version 1.1 fixes full-sample-ID detection, including mixed IDs such as
`170726001-Hg-59` and `090726D01-Hg-44`. Repeated readings use the full ID, so
the trailing number is not mistaken for a replicate and non-D-token IDs are
not discarded. The regression, ESD, and CI calculations are unchanged.
Small numerical fits use one computational thread to avoid unnecessary thread
overhead. This changes execution speed, not the statistical method.

## Run

Use Python 3.12 (the supplied requirements retain the original app's pinned versions).
Open a terminal in this folder:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Windows: double-click `run_local.bat` after installing Python. macOS/Linux:
`bash run_local.sh`. The launchers install the requirements if needed.

For Streamlit Community Cloud, upload **all four** of `app.py`,
`regression_core.py`, `retained_methods.py`, and `ui_controls.py`, plus `requirements.txt` and
`sample_input_template.csv`. Select `app.py` as the main file. `run_app.py` and
`streamlit_app.py` are alternative launch entry points. This package has not been
deployed to an external hosting service in the original ZIP.

### Streamlit Community Cloud deployment

Use repository `Abicumaran/Regression-2`, branch `main`, and main file `app.py`.
Before deploying, open **Advanced settings**, select **Python 3.12**, and save.
The pinned NumPy, pandas, and SciPy releases in this project are not intended for
Python 3.14. A newer interpreter can trigger failed source builds instead of
installing compatible prebuilt wheels. Installing a compiler is not the fix.

The `.python-version` file records the local runtime target; it does **not**
replace the Python version selection in Streamlit Community Cloud. For an
existing deployment on the wrong Python version, follow Streamlit's documented
delete-and-redeploy procedure, preserving the repository and deployment settings:
https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/upgrade-python

## Input

Upload `.xlsx` or `.csv`; select the worksheet for Excel workbooks. The first row
must contain unique column headings, followed by one paired replicate per row.
Save legacy `.xls` files as `.xlsx` first.

- Measured/reference pairs such as `RBC` / `RBC_ref`, `PLT_2` / `PLT_ref`, and
  `PLT_3` / `PLT_ref` are detected automatically. Other markers with `_ref`,
  `_reference`, or recognized reference prefixes/suffixes work too. Use the
  column-mapping selectors for arbitrary X and Y names. Each model has its own
  result and can share a genuine reference column.
- Donor, specimen, replicate, sample ID, batch ID, and global flag are detected
  where available. Review the mappings in the app. Original sample IDs such as
  `260101D01-Mc-1` are supported. Explicit donor labels take priority. A consistent
  recognized donor/specimen/replicate layout uses collection prefix + D-token.
  Other or mixed ID layouts use the full sample ID unchanged. Numeric suffixes
  conflicting with the explicit replicate column are not stripped automatically.
- If no donor or sample-ID column is present, each row is an independent sample.
  If those rows are repeated measurements, select an appropriate donor column.
- If a sample-ID column is present but lacks a recognized replicate layout,
  identical full IDs are averaged as repeat readings. Confidence intervals then
  assume those sample IDs are independent; supply a donor column if different
  sample IDs belong to one donor. `analysis_unit` and `n_independent_units` identify
  the grouping used. `n_donors` and the `Donor Means` sheet retain their legacy names.
- If no specimen type is available, unstructured inputs use one `All` group.
  Structured sample IDs retain parsed specimen groups. Mc/Mc2 are combined as in
  the original app. Groups and donor selections are editable.
- If no global-flag column is selected, rows are treated as FALSE and the app
  states this. With a flag column, TRUE and missing/unrecognized flags are excluded.
- Use measured and reference values in compatible units. Only complete finite
  pairs are retained for each analyte. No reference values are fabricated.
- Optional normal-range plot markers retain the original defaults for known
  analytes. Change or clear those defaults for another population. These markers
  never filter observations and are separate from normality testing.

`sample_input_template.csv` contains **synthetic demonstration data**: 24 donors,
two specimen types, three replicates per donor/type, and two measured/reference
pairs. It includes a flagged row, an extreme value, and a missing value for
demonstrating the audits. It is not study data.

## Statistical methods

The original app had different outlier logic in its two analysis branches.
Regression 2 retains the **regression branch**: global flag FALSE → original
generalized ESD on raw linked replicate rows → analyte-specific removal → paired
donor means → Huber regression. It contains no Bland–Altman analysis or report.

The 19 retained functions in `retained_methods.py` are copied unchanged from the
source app's **embedded** engine (`2026-09-27-proxima-v7-validation`), which is newer
than the standalone `analysis_core.py` in that ZIP.

| Quantity | Method |
| --- | --- |
| Slope and intercept | Original `HuberRegressor()` defaults; X = reference, Y = measured |
| Slope/intercept 95% CIs | Original paired bootstrap of donor means; 2.5th/97.5th percentiles; seed 20260722 |
| Pearson R and p | Original Pearson correlation on the same retained paired donor means |
| R 95% CI | Original Fisher-z interval; SE = 1/sqrt(n − 3) |
| Outliers | Original generalized ESD with SD/CV/mixed error-pattern routing and original +10000 log transformation |
| Normality | Shapiro–Wilk diagnostics at alpha 0.05 for raw paired differences and final Huber residuals |

All three 95% CIs are requested automatically. The original app's optional bootstrap
toggle is removed; the default is now 5,000 resamples, adjustable from 200 to 20,000.
The bootstrap resamples paired donor means, **not individual replicates**. ESD and
averaging are fixed before resampling, as in the source, so the intervals do not
include uncertainty in outlier selection. Successful bootstrap counts and fitting
warnings are reported.

Normality diagnostics and residual/Q-Q plots do not switch the regression estimator
or ESD method. The original Shapiro → Grubbs/MAD removal rule belonged to the
separate agreement branch, so it is not substituted for regression ESD here.

The output also distinguishes:

- Pearson R squared, which is the squared correlation.
- R squared and RMSE of the fitted Huber line, comparing observed Y to fitted Y.
- The source app's R squared and RMSE against identity, comparing Y directly to X.

At least four independent donor pairs are required for the requested CIs. Constant
X/Y or insufficient data produce unavailable results with an explicit status.
Generalized ESD is not applied to tested subsets with fewer than 20 complete rows,
matching the original implementation. The inherited transformation requires
values greater than −10001 when ESD is needed; an incompatible analyte is explicitly
marked not analyzed. Finite values, unique internal row IDs, and separate pair
processing prevent missing values, duplicate batch labels, or shared references
from silently changing the wrong records.

## Downloads

- `Regression_2_Results.xlsx`: regression estimates and all CI bounds, donor means,
  outlier and exclusion audits, ESD diagnostics, normality, mappings, and methods.
- `Regression_2_Metrics.csv`: the full regression metrics table.
- `Regression_2_Plots.zip`: individual regression and diagnostic plots, plus
  combined regression panels, as 300-dpi PNG and vector PDF files.

Plots are generated on request after the analysis. Changed input/settings hide
stale results until the analysis is rerun. The app writes no uploaded data to an
external service itself; data are processed in the running Streamlit server.

## Verification

Run `python -m unittest discover -s tests -v` from this directory.
See `VALIDATION.txt` for the completed source-parity and workflow checks.

Reference documentation:

- https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.HuberRegressor.html
- https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats._result_classes.PearsonRResult.confidence_interval.html
