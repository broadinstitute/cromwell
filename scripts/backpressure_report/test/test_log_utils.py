import unittest
from backpressure_report.lib import log_utils
from test.lib.log_helper import LOG


class LogUtilsTestMethods(unittest.TestCase):

    def test_build_log_jsons_from_input_files(self):
        self.assertEqual(34, len(LOG))

    def test_is_start_event(self):
        starts = [e for e in LOG if log_utils.is_start_event(e)]
        self.assertEqual(17, len(starts))

    def test_is_end_event(self):
        ends = [e for e in LOG if log_utils.is_end_event(e)]
        self.assertEqual(17, len(ends))

    def test_filter_and_sort_log_entries(self):
        filtered = log_utils.filter_and_sort_log_entries(LOG)
        self.assertEqual(34, len(filtered))
        oldest = filtered[0]['timestamp']
        self.assertTrue(all([t['timestamp'] >= oldest for t in filtered]))
        self.assertTrue(all(log_utils.is_start_event(e) or log_utils.is_end_event(e) for e in filtered))
