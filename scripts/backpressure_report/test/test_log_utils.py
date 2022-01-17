import pathlib
import unittest
from backpressure_report.lib import log_utils

__log_path = pathlib.Path(__file__).parent / 'resources' / 'alpha_logs.json'
LOG = log_utils.build_log_jsons_from_input_files([__log_path])[0]


class LogUtilsTestMethods(unittest.TestCase):

    def test_build_log_jsons_from_input_files(self):
        self.assertEqual(34, len(LOG))

    def test_is_event_start(self):
        starts = [e for e in LOG if log_utils.is_event_start(e)]
        self.assertEqual(17, len(starts))

    def test_is_event_end(self):
        starts = [e for e in LOG if log_utils.is_event_end(e)]
        self.assertEqual(17, len(starts))

    def test_filter_and_sort_log_entries(self):
        filtered = log_utils.filter_and_sort_log_entries(LOG)
        self.assertEqual(34, len(filtered))
        oldest = filtered[0]['timestamp']
        self.assertTrue(all([t['timestamp'] >= oldest for t in filtered]))
