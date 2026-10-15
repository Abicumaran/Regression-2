"""Bias formulas, retained populations, and report exports."""
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import regression_core as core


class BiasTests(unittest.TestCase):
    def test_screenshot_formula_and_signed_direction(self):
        result = core.bias_metrics([10, 20, 30], [12, 15, 27])
        self.assertEqual(result["mean_true"], 20)
        self.assertEqual(result["bias"], -2)
        self.assertEqual(result["bias_percent"], -10)
        self.assertEqual(result["abs_bias_percent"], 10)
        self.assertEqual(result["bias_status"], "OK")
        # Opposite signed errors cancel; this is not mean absolute error.
        self.assertEqual(core.bias_metrics([10, 10], [8, 12])["abs_bias_percent"], 0)

    def test_pairwise_finite_zero_mean_and_empty(self):
        result = core.bias_metrics([10, np.nan, 20, np.inf], [12, 100, np.inf, 8])
        self.assertEqual(result["n_bias_pairs"], 1)
        self.assertEqual(result["bias_percent"], 20)
        result = core.bias_metrics([-1, 1], [0, 2])
        self.assertEqual(result["bias"], 1)
        self.assertTrue(np.isnan(result["bias_percent"]))
        self.assertIn("zero", result["bias_status"])
        result = core.bias_metrics([], [])
        self.assertEqual(result["n_bias_pairs"], 0)
        self.assertTrue(np.isnan(result["bias"]))
        with self.assertRaises(ValueError):
            core.bias_metrics([1, 2], [1])

    def analyze(self, data):
        prepared, _ = core.prepare_metadata(data, core.detect_metadata(data))
        return core.analyze(data, prepared, core.detect_pairs(data.columns), {"All": ["All"]}, n_boot=5)

    def test_retained_readings_and_donor_means_use_different_weights(self):
        data = pd.DataFrame({"Donor": ["A", "A", "B", "C", "D", "E"],
                             "Marker_ref": [10, 10, 20, 100, 30, 40],
                             "Marker": [12, 12, 19, 200, np.nan, 99],
                             "global_flag": [False, False, False, True, False, False]})
        def remove_last(work, group, analytes):
            return pd.DataFrame({"analyte": ["Marker"], "batch_id": [work.iloc[-1].batch_id]}), pd.DataFrame()
        with patch.object(core.original, "detect_esd_outliers", side_effect=remove_last):
            result = self.analyze(data)
        row = result["tables"]["Regression Results"].iloc[0]
        self.assertEqual(row.n_bias_pairs, 3)
        self.assertEqual(row.bias, 1)
        self.assertAlmostEqual(row.bias_percent, 7.5)
        self.assertEqual(row.donor_mean_n_bias_pairs, 2)
        self.assertEqual(row.donor_mean_bias, .5)
        self.assertAlmostEqual(row.donor_mean_bias_percent, 100 * .5 / 15)
        sheets = pd.read_excel(io.BytesIO(core.export_workbook(result)), sheet_name=None)
        self.assertIn("Bias Results", sheets)
        for name in ["Bias Results", "Regression Results"]:
            self.assertEqual(sheets[name].iloc[0].bias, 1)
            self.assertAlmostEqual(sheets[name].iloc[0].abs_bias_percent, 7.5)
        csv = pd.read_csv(io.StringIO(result["tables"]["Regression Results"].to_csv(index=False)))
        self.assertEqual(csv.iloc[0].bias, 1)

    def test_bias_available_without_regression_fit(self):
        for refs, measured in [([10], [12]), ([10, 10], [12, 12])]:
            row = self.analyze(pd.DataFrame({"Marker_ref": refs, "Marker": measured}))["tables"]["Regression Results"].iloc[0]
            self.assertEqual(row.bias, 2)
            self.assertEqual(row.bias_percent, 20)
            self.assertEqual(row.bias_status, "OK")
            self.assertTrue(np.isnan(row.huber_slope))


if __name__ == "__main__":
    unittest.main()
