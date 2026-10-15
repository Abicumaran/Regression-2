"""Exercise the deployed entry point and synthetic-data workflow on Python 3.12."""
import io
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch
import sys

import pandas as pd
from streamlit.testing.v1 import AppTest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ui_controls import ANALYTE_ORDER, analyte_choices


class AppWorkflowTests(unittest.TestCase):
    def test_requested_order_missing_columns_and_extra_analytes(self):
        choices = analyte_choices(["Extra_ref", "MXD_ref", "PLT_ref", "MXD_2", "Extra", "PLT_3"])
        self.assertEqual([label for label, _ in choices[:12]], list(ANALYTE_ORDER))
        self.assertEqual(choices[-1][0], "Extra")
        self.assertFalse(choices[0][1]["Run"])
        self.assertIsNone(choices[0][1]["Measured (Y)"])
        self.assertTrue(choices[2][1]["Run"])
        self.assertEqual(choices[2][1]["Reference (X)"], "PLT_ref")

    def test_order_exclusion_and_plain_tables_in_app(self):
        data = {"Donor": [f"D{i}" for i in range(8)]}
        for name in reversed(ANALYTE_ORDER):
            data[name] = [1.2 * i + 0.4 + (i % 2) * 0.01 for i in range(1, 9)]
            data[name + "_ref"] = list(range(1, 9))
        frame = pd.DataFrame(data)
        with patch("regression_core.read_upload", return_value={"CSV": frame}):
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=90).run()
            app.checkbox[0].check().run()
            analytes = [c for c in app.checkbox if c.label.startswith("Include ") and not c.label.startswith("Include group:")]
            self.assertEqual([c.label for c in analytes], [f"Include {name}" for name in ANALYTE_ORDER])
            self.assertTrue(all(c.value for c in analytes))
            for checkbox in analytes:
                if checkbox.label not in ("Include RBC_3", "Include PLT_3", "Include MXD_2"):
                    checkbox.uncheck()
            next(n for n in app.number_input if n.label == "Donor-level bootstrap resamples").set_value(200)
            next(b for b in app.button if b.label == "Run Regression 2").click().run()
            self.assertEqual(len(app.exception), 0)
            self.assertEqual(len(app.error), 0)
            metrics = app.session_state["regression_2_output"]["result"]["tables"]["Regression Results"]
            self.assertEqual(metrics.analyte.tolist(), ["RBC_3", "PLT_3", "MXD_2"])
            self.assertEqual(len(app.get("arrow_data_frame")), 0)
            self.assertTrue(any('r2-table-wrap' in element.value for element in app.markdown))
            next(c for c in app.checkbox if c.label == "Include PLT_3").uncheck().run()
            self.assertTrue(any("Settings or input changed" in item.value for item in app.info))

    def test_manual_mapping_missing_analyte_and_custom_group(self):
        frame = pd.DataFrame({"Donor": ["A", "B", "C", "D"], "X": [1, 2, 3, 4], "Y": [2.1, 4.2, 6.1, 8.3]})
        with patch("regression_core.read_upload", return_value={"CSV": frame}):
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=90).run()
            app.checkbox[0].check().run()
            next(c for c in app.checkbox if c.label == "Include RBC_3").check()
            next(s for s in app.selectbox if s.label == "Measured (Y)").select("Y")
            next(s for s in app.selectbox if s.label == "Reference (X)").select("X")
            next(n for n in app.number_input if n.label == "Additional custom groups").set_value(1).run()
            next(c for c in app.checkbox if c.label == "Include group: All").uncheck()
            next(c for c in app.checkbox if c.label == "Include group: Custom group 1").check()
            [s for s in app.multiselect if s.label == "Specimen types"][-1].set_value(["All"])
            next(n for n in app.number_input if n.label == "Donor-level bootstrap resamples").set_value(200)
            next(b for b in app.button if b.label == "Run Regression 2").click().run()
            self.assertEqual(len(app.exception), 0)
            self.assertEqual(len(app.error), 0)
            metrics = app.session_state["regression_2_output"]["result"]["tables"]["Regression Results"]
            self.assertEqual(metrics.analyte.tolist(), ["Y"])
            self.assertEqual(metrics.analysis_group.tolist(), ["Custom group 1"])

    def test_synthetic_analysis_and_downloads(self):
        app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"),
                               default_timeout=90).run()
        self.assertEqual(len(app.exception), 0)
        app.checkbox[0].check().run()
        self.assertEqual(len(app.exception), 0)
        next(n for n in app.number_input if n.label == "Donor-level bootstrap resamples").set_value(200)
        next(b for b in app.button if b.label == "Run Regression 2").click().run()
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(len(app.error), 0)
        saved = app.session_state["regression_2_output"]
        metrics = saved["result"]["tables"]["Regression Results"]
        self.assertFalse(metrics.empty)
        self.assertTrue(metrics["huber_slope"].notna().all())
        workbook = pd.read_excel(io.BytesIO(saved["workbook"]), sheet_name=None)
        self.assertIn("Regression Results", workbook)
        self.assertIn("Bias Results", workbook)
        self.assertTrue(metrics.bias.notna().all())
        self.assertTrue(any("Bias — retained paired readings" in item.value for item in app.subheader))
        self.assertTrue(any("Absolute bias (%)" in item.value for item in app.markdown))
        next(b for b in app.button if b.label.startswith("Prepare all plots")).click().run()
        self.assertEqual(len(app.exception), 0)
        with zipfile.ZipFile(io.BytesIO(app.session_state["regression_2_plots"])) as plots:
            self.assertTrue(any(name.endswith(".png") for name in plots.namelist()))
            self.assertTrue(any(name.endswith(".pdf") for name in plots.namelist()))
        next(n for n in app.number_input if n.label == "Donor-level bootstrap resamples").set_value(400).run()
        self.assertTrue(any("Settings or input changed" in item.value for item in app.info))


if __name__ == "__main__":
    unittest.main()
