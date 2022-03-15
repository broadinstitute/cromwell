import unittest

from backpressure_report.lib import backpressure_event, backpressure_window
from test.lib.log_helper import LOG


class BackpressureWindowTestMethods(unittest.TestCase):

    def setUp(self) -> None:
        self.events = backpressure_event.build_backpressure_events_from_log_jsons([LOG])
        self.expected_durations = [36, 0, 32, 0, 189, 0, 0, 0, 46, 0, 0, 142, 60, 29]

    def __run_assert(self) -> None:
        windows, all_pods = backpressure_window.build_windows_and_pods_from_events(self.events, window_width_in_hours=8)
        actual_durations = [
            sum(w.durations_by_pod().values()) for w in windows
        ]
        self.assertEqual(self.expected_durations, actual_durations)
        # One pod backpressuring
        self.assertEqual("2022-01-11 00:00:00+00:00,32,0,0,0,0,0,32,0,0", windows[2].report_line(all_pods))
        # No pods backpressuring
        self.assertEqual("2022-01-11 08:00:00+00:00,0,0,0,0,0,0,0,0,0", windows[3].report_line(all_pods))
        # Multiple pods backpressuring
        self.assertEqual("2022-01-11 16:00:00+00:00,189,76,0,89,0,0,0,24,0", windows[4].report_line(all_pods))

    def test_forward(self):
        self.__run_assert()

    def test_reversed(self):
        self.events.reverse()
        self.__run_assert()
