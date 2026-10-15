"""Numerical and input-safety checks. Run: python -m unittest discover -s tests -v."""
import io
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import regression_core as core
import retained_methods as original


class RegressionTests(unittest.TestCase):
    def data(self, n=12, reps=2):
        rng = np.random.default_rng(17)
        x = np.repeat(np.linspace(1, 8, n), reps)
        return pd.DataFrame({"Donor": np.repeat([f"D{i:02}" for i in range(n)], reps),
            "global_flag": False, "Marker": 1.2*x+.4+rng.normal(0,.08,len(x)), "Marker_ref": x})

    def run_data(self, data, boot=20):
        prep, _ = core.prepare_metadata(data, core.detect_metadata(data))
        pairs = core.validate_pairs(core.detect_pairs(data.columns), data.columns)
        return core.analyze(data, prep, pairs, {"All": ["All"]}, n_boot=boot)

    def test_fisher_matches_scipy(self):
        x = np.arange(1, 21, dtype=float)
        y = 1.3*x+np.random.default_rng(5).normal(size=20)
        fit = stats.pearsonr(x, y)
        expected = fit.confidence_interval(confidence_level=.95)
        np.testing.assert_allclose(original.pearson_ci(fit.statistic, len(x)), [expected.low, expected.high], atol=1e-14)

    def test_known_line_and_fit_metrics(self):
        x = np.arange(1, 21, dtype=float)
        row = core.regression_metrics(x, 2*x+3, "All", "Y", "X", "", 20)
        self.assertAlmostEqual(row["huber_slope"], 2, places=4)
        self.assertAlmostEqual(row["huber_intercept"], 3, places=4)
        self.assertAlmostEqual(row["huber_fit_r_squared"], 1, places=7)
        self.assertLess(row["huber_fit_rmse"], 1e-4)
        self.assertGreater(row["rmse_vs_identity"], 10)

    def test_averaging_and_flags_and_missing_are_pairwise(self):
        d = self.data(n=8)
        d["global_flag"] = d["global_flag"].astype(object)
        d.loc[0, "global_flag"] = True
        d.loc[2, "global_flag"] = None
        d.loc[4, "Marker"] = np.nan
        result = self.run_data(d)
        means = result["tables"]["Donor Means"].set_index("Donor")
        self.assertEqual(len(means), 8)
        for donor in ["D00", "D01", "D02"]:
            self.assertEqual(means.loc[donor, "retained_replicates"], 1)
        self.assertEqual(len(result["tables"]["global flag exclusions"]), 2)
        self.assertEqual(len(result["tables"]["Missing Pairs"]), 1)
        self.assertEqual(result["tables"]["Regression Results"].iloc[0].n_donors, 8)

    def test_shared_reference_model_variants(self):
        d = self.data(n=8).rename(columns={"Marker": "Marker_2"})
        d["Marker_3"] = d["Marker_2"] + 1
        result = self.run_data(d)
        m = result["tables"]["Regression Results"].set_index("analyte")
        self.assertEqual(set(m.index), {"Marker_2", "Marker_3"})
        self.assertEqual(set(m.reference_column), {"Marker_ref"})
        self.assertAlmostEqual(m.loc["Marker_3","huber_intercept"]-m.loc["Marker_2","huber_intercept"], 1, places=3)

    def test_small_and_constant_data_have_explicit_status(self):
        d = self.data(n=3)
        m = self.run_data(d)["tables"]["Regression Results"].iloc[0]
        self.assertTrue(np.isnan(m.huber_slope_ci95_low))
        self.assertTrue(np.isnan(m.pearson_r_ci95_low))
        self.assertIn("at least 4", m.status)
        d["Marker_ref"] = 1
        m = self.run_data(d)["tables"]["Regression Results"].iloc[0]
        self.assertIn("constant", m.status)

    def test_identity_mapping_and_unique_rows(self):
        d = pd.DataFrame({"bloodSampleId": ["260101D01-Mc-1", "260201D01-Mc-1"],
                          "batch_id": ["duplicate", "duplicate"]})
        prep, notes = core.prepare_metadata(d, core.detect_metadata(d))
        self.assertEqual(prep.donor.nunique(), 2)
        self.assertEqual(prep.batch_id.nunique(), 2)
        self.assertTrue(prep.global_flag_bool.eq(False).all())
        with self.assertRaises(ValueError):
            core.prepare_metadata(d, core.detect_metadata(d), "token")

    def test_custom_names_and_no_metadata(self):
        d = pd.DataFrame({"X": [1,2,3,4,5], "Y": [2.1,3.8,6.2,7.9,10.3]})
        mapping = core.detect_metadata(d)
        prep, notes = core.prepare_metadata(d, mapping)
        self.assertEqual(prep.donor.nunique(), 5)
        pairs = pd.DataFrame([{"Run": True, "Measured (Y)": "Y", "Reference (X)": "X", "Unit": ""}])
        pairs = core.validate_pairs(pairs, d.columns)
        result = core.analyze(d, prep, pairs, {"All": ["All"]}, n_boot=20)
        self.assertEqual(result["tables"]["Regression Results"].iloc[0].n_donors, 5)
        pairs.loc[0, "Reference (X)"] = "Y"
        with self.assertRaises(ValueError):
            core.validate_pairs(pairs, d.columns)

    def test_generic_sample_ids_with_explicit_donor(self):
        d = self.data(n=8)
        d["Sample ID"] = [f"tube-{i}" for i in range(len(d))]
        prep, _ = core.prepare_metadata(d, core.detect_metadata(d))
        self.assertTrue(prep.parse_ok.all())
        self.assertEqual(prep.donor.nunique(), 8)

    def test_bootstrap_repeatability_and_bounds(self):
        d = self.data(n=8)
        first = self.run_data(d)["tables"]["Regression Results"]
        second = self.run_data(d)["tables"]["Regression Results"]
        pd.testing.assert_frame_equal(first, second)
        m = first.iloc[0]
        for kind in ["huber_slope", "huber_intercept", "pearson_r"]:
            self.assertLessEqual(m[kind+"_ci95_low"], m[kind+"_ci95_high"])
        self.assertEqual(m.bootstrap_successful, m.bootstrap_requested)

    def test_transform_domain_is_explicit(self):
        d = self.data(n=12)
        d.loc[0, "Marker"] = -20000
        m = self.run_data(d)["tables"]["Regression Results"].iloc[0]
        self.assertIn("requires values > -10001", m.status)

    def test_excel_export_and_excel_input(self):
        result = self.run_data(self.data(n=8))
        payload = core.export_workbook(result)
        sheets = pd.read_excel(io.BytesIO(payload), sheet_name=None)
        self.assertIn("Regression Results", sheets)
        self.assertIn("Normality", sheets)
        self.assertFalse(any("bland" in s.lower() for s in sheets))
        m = sheets["Regression Results"]
        self.assertTrue(m.huber_slope_ci95_low.notna().all())
        loaded = core.read_upload("report.xlsx", payload)
        self.assertEqual(len(loaded["Regression Results"]), len(m))

    def test_group_with_all_rows_flagged_does_not_stop_other_groups(self):
        d = self.data(n=8)
        d["Type"] = "A"
        extra = d.iloc[[0]].copy()
        extra["Type"], extra["global_flag"] = "B", True
        d = pd.concat([d, extra], ignore_index=True)
        prep, _ = core.prepare_metadata(d, core.detect_metadata(d))
        result = core.analyze(d, prep, core.detect_pairs(d.columns), {"A": ["A"], "B": ["B"]}, n_boot=20)
        metrics = result["tables"]["Regression Results"].set_index("analysis_group")
        self.assertEqual(metrics.loc["A", "n_donors"], 8)
        self.assertEqual(metrics.loc["B", "n_donors"], 0)
        self.assertIn("Unavailable", metrics.loc["B", "status"])

    def test_mixed_hg_ids_keep_full_identity_and_all_rows(self):
        ids = ["170726001-Hg-59", "170726001-Hg-140", "090726D01-Hg-44", "090726D01-Hg-142"]
        d = pd.DataFrame({"bloodSampleId": np.repeat(ids, 2), "Replicate": 1,
                          "global_flag": False, "X": np.repeat([1,2,3,4], 2), "Y": np.arange(8)+1.0})
        prep, _ = core.prepare_metadata(d, core.detect_metadata(d))
        self.assertTrue(prep.parse_ok.all())
        self.assertEqual(set(prep.donor), set(ids))
        self.assertTrue(prep.donor_identity_rule.eq("full_sample").all())
        self.assertTrue(prep.specimen_type.eq("All").all())
        self.assertEqual(prep.bloodSampleId.tolist(), d.bloodSampleId.tolist())

    def test_unknown_tag_numbers_conflicting_with_replicate_are_not_stripped(self):
        d = pd.DataFrame({"bloodSampleId": ["090726D01-Hg-44", "090726D01-Hg-142"], "Replicate": [1,1]})
        prep, _ = core.prepare_metadata(d, core.detect_metadata(d))
        self.assertEqual(prep.donor.nunique(), 2)
        self.assertTrue(prep.donor_identity_rule.eq("full_sample").all())

    def test_original_structured_ids_still_average_replicates(self):
        d = pd.DataFrame({"bloodSampleId": ["260101D01-Mc-1", "260101D01-Mc-2", "260101D02-Mc-1"], "Replicate": [1,2,1]})
        prep, _ = core.prepare_metadata(d, core.detect_metadata(d))
        self.assertEqual(prep.donor.nunique(), 2)
        self.assertTrue(prep.donor_identity_rule.eq("prefix").all())
        self.assertTrue(prep.specimen_type.eq("Mc").all())


if __name__ == "__main__":
    unittest.main()
