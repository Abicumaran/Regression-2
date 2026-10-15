"""Exercise interrupted page reruns and failures without losing completed results."""
import unittest
from pathlib import Path
from threading import Event
from unittest.mock import patch
import sys

import pandas as pd
from streamlit.testing.v1 import AppTest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import regression_core
from test_app import finish_job


class PersistenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = pd.read_csv(Path(__file__).resolve().parents[1] / "sample_input_template.csv")
        prepared, _ = regression_core.prepare_metadata(data, regression_core.detect_metadata(data), "auto")
        pairs = regression_core.validate_pairs(regression_core.detect_pairs(data.columns), data.columns)
        groups = regression_core.validate_groups(regression_core.default_groups(sorted(prepared.specimen_type.unique())), sorted(prepared.specimen_type.unique()))
        cls.result = regression_core.analyze(data, prepared, pairs, groups, n_boot=200)

    def app(self):
        app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=30).run()
        app.checkbox[0].check().run()
        next(n for n in app.number_input if n.label == "Donor-level bootstrap resamples").set_value(200)
        return app

    def start(self, app):
        next(b for b in app.button if b.label == "Run Regression 2").click().run()

    def assert_results(self, app):
        self.assertEqual(len(app.exception), 0)
        self.assertTrue(any(x.value == "Regression results" for x in app.subheader))
        self.assertTrue(any("Analysis complete" in x.value for x in app.success))
        self.assertTrue(any(x.label == "Download all regression metrics (CSV)" for x in app.get("download_button")))

    def test_running_job_survives_control_changes_and_input_removal(self):
        release = Event()
        received = {}

        def delayed(*args, **kwargs):
            received["n_boot"] = args[5]
            kwargs["progress"](0.5, "Testing progress")
            if not release.wait(20):
                raise RuntimeError("Test release timed out")
            return {**self.result, "tables": {k: v.copy() for k, v in self.result["tables"].items()}}

        with patch("analysis_job.analyze", side_effect=delayed) as analyze:
            app = self.app()
            try:
                self.start(app)
                job = app.session_state["regression_2_job"]
                self.assertTrue(next(b for b in app.button if b.label == "Run Regression 2").disabled)
                next(n for n in app.number_input if n.label == "Donor-level bootstrap resamples").set_value(400).run()
                self.assertIs(app.session_state["regression_2_job"], job)
                app.checkbox[0].uncheck().run()
                self.assertTrue(any("Analysis is running" in x.value for x in app.info))
            finally:
                release.set()
            finish_job(app)
            self.assertEqual(analyze.call_count, 1)
            self.assertEqual(received["n_boot"], 200)
            self.assert_results(app)
            self.assertIsNotNone(app.session_state["regression_2_output"]["workbook"])

    def test_export_failure_and_failed_rerun_keep_results(self):
        with patch("analysis_job.analyze", return_value=self.result), patch("analysis_job.export_workbook", side_effect=RuntimeError("test export failure")):
            app = self.app()
            self.start(app)
            finish_job(app)
        self.assert_results(app)
        self.assertIsNone(app.session_state["regression_2_output"]["workbook"])
        self.assertTrue(any("Excel report could not" in x.value for x in app.warning))
        next(b for b in app.button if b.label == "Retry Excel report").click().run()
        self.assertIsNotNone(app.session_state["regression_2_output"]["workbook"])
        prior = app.session_state["regression_2_output"]
        with patch("analysis_job.analyze", side_effect=RuntimeError("test analysis failure")):
            self.start(app)
            finish_job(app)
        self.assert_results(app)
        self.assertIs(app.session_state["regression_2_output"], prior)
        self.assertTrue(any("test analysis failure" in x.value for x in app.error))

    def test_plot_error_does_not_hide_report_or_tables(self):
        with patch("analysis_job.analyze", return_value=self.result):
            app = self.app()
            self.start(app)
            finish_job(app)
        with patch("results_ui.regression_figure", side_effect=RuntimeError("test plot failure")):
            app.run()
        self.assert_results(app)
        self.assertTrue(any("Plot could not be prepared" in x.value for x in app.error))
        self.assertTrue(any(x.label == "Download Regression 2 Excel report" for x in app.get("download_button")))


if __name__ == "__main__":
    unittest.main()
