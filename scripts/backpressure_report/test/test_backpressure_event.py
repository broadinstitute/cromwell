import unittest

from backpressure_report.lib import backpressure_event
from test.lib.log_helper import LOG


class BackpressureEventTestMethods(unittest.TestCase):

    def test_build_backpressure_events_from_log_jsons(self):
        events = backpressure_event.build_backpressure_events_from_log_jsons([LOG])
        self.assertEqual(17, len(events))
        # Minimum backpressure duration is 20 seconds, maximum is 60 seconds.
        self.assertTrue(all([(d.duration() >= 20) and (d.duration() <= 60) for d in events]))
