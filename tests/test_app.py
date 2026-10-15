"""Exercise the deployed entry point and synthetic-data workflow on Python 3.12."""
import io
import unittest
import zipfile
from pathlib import Path

import pandas as pd
from streamlit.testing.v1 import AppTest


class AppWorkflowTests(unittest.TestCase):
    def test_synthetic_analysis_and_downloads(self):
        app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"),
                               default_timeout=90).run()
        self.assertEqual(len(app.exception), 0)
        app.checkbox[0].check().run()
        self.assertEqual(len(app.exception), 0)
        app.number_input[0].set_value(200)
        next(b for b in app.button if b.label == "Run Regression 2").click().run()
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(len(app.error), 0)
        saved = app.session_state["regression_2_output"]
        metrics = saved["result"]["tables"]["Regression Results"]
        self.assertFalse(metrics.empty)
        self.assertTrue(metrics["huber_slope"].notna().all())
        workbook = pd.read_excel(io.BytesIO(saved["workbook"]), sheet_name=None)
        self.assertIn("Regression Results", workbook)
        next(b for b in app.button if b.label.startswith("Prepare all plots")).click().run()
        self.assertEqual(len(app.exception), 0)
        with zipfile.ZipFile(io.BytesIO(app.session_state["regression_2_plots"])) as plots:
            self.assertTrue(any(name.endswith(".png") for name in plots.namelist()))
            self.assertTrue(any(name.endswith(".pdf") for name in plots.namelist()))
        app.number_input[0].set_value(400).run()
        self.assertTrue(any("Settings or input changed" in item.value for item in app.info))


if __name__ == "__main__":
    unittest.main()
