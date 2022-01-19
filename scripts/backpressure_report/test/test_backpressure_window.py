import unittest

from backpressure_report.lib import backpressure_event, backpressure_window
from test.lib.log_helper import LOG


class BackpressureWindowTestMethods(unittest.TestCase):

    def setUp(self) -> None:
        self.events = backpressure_event.build_backpressure_events_from_log_jsons([LOG])
        self.expected = [36, 0, 32, 0, 189, 0, 0, 0, 46, 0, 0, 142, 60, 29]

    def __actual(self) -> list:
        windows, all_pods = backpressure_window.build_windows_from_events(self.events, window_width_in_hours=8)
        return [
            sum(w.durations_by_pod().values()) for w in windows
        ]

    def test_forward(self):
        self.assertEqual(self.expected, self.__actual())

    def test_reversed(self):
        self.events.reverse()
        self.assertEqual(self.expected, self.__actual())
