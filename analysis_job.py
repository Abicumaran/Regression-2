"""Session-owned analysis workers; no Streamlit calls or shared result cache."""
from datetime import datetime, timezone
from threading import Lock, Thread

from regression_core import analyze, export_workbook


class AnalysisJob:
    def __init__(self, data, prepared, pairs, groups, donors, n_boot, *, signature,
                 source, show_ranges):
        self._lock = Lock()
        self._state = dict(done=False, progress=0.0, message="Starting regression analysis…",
                           output=None, error=None)
        self.handled = False
        self.thread = Thread(target=self._run, args=(data.copy(deep=True), prepared.copy(deep=True),
            pairs.copy(deep=True), {key: list(value) for key, value in groups.items()},
            list(donors), n_boot, signature, source, show_ranges), daemon=True)

    def start(self):
        self.thread.start()

    def snapshot(self):
        with self._lock:
            return self._state.copy()

    def _update(self, **values):
        with self._lock:
            self._state.update(values)

    def _progress(self, value, message):
        self._update(progress=float(value), message=message)

    def _run(self, data, prepared, pairs, groups, donors, n_boot, signature, source, show_ranges):
        try:
            result = analyze(data, prepared, pairs, groups, donors, n_boot, progress=self._progress)
            completed = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
            methods = result["tables"]["Methods"]
            methods.loc[len(methods)] = ["Input", source]
            methods.loc[len(methods)] = ["Completed", completed]
            output = dict(signature=signature, result=result, workbook=None,
                          completed=completed, source=source, show_ranges=show_ranges,
                          export_error=None)
            # Publish the completed analysis before attempting an optional file export.
            self._update(output=output, progress=1.0, message="Analysis complete. Preparing Excel report…")
            try:
                output = {**output, "workbook": export_workbook(result)}
            except Exception as exc:
                output = {**output, "export_error": str(exc)}
            self._update(output=output)
        except Exception as exc:
            self._update(error=str(exc))
        finally:
            self._update(done=True)
