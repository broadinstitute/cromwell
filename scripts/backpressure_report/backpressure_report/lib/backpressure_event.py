from backpressure_report.lib import log_utils
from datetime import datetime
from dateutil import parser
import itertools
from typing import AnyStr


class BackpressureEvent:
    """
    Represents a span during which a single pod was backpressuring / in high I/O and not dispensing job tokens.
    """
    def __init__(self, pod: str, start: datetime, end: datetime):
        self.pod = pod
        self.start = start
        self.end = end

    def duration(self):
        return (self.end - self.start).seconds

    def __str__(self) -> AnyStr:
        return f"BackpressureEvent(pod = {self.pod},start={str(self.start)},duration={(self.duration())}s)"


def build_backpressure_events_from_log_jsons(logs: list):
    """
    Build a list of BackpressureEvents from the specified logs, using matched "start" and "end" events for a particular
    pod to delimit the duration of the BackpressureEvent.

    :param logs: a list of JSON log files, each of which is a list of JSON objects each representing a log entry.
    :return: a list of BackpressureEvents.
    """

    # Complete BackpressureEvent objects corresponding to a matched pair of backpressure start and stop log entries for
    # a pod.
    complete = []
    # Already-processed log entry ids to ignore duplicates in overlapping log file ranges.
    seen_insert_ids = set(())
    # pod names for which we have seen a "backpressure start" log messages and for which we are now awaiting a matching
    # "backpressure stop" log message for the same pod name.
    in_progress_starts_by_pod = {}

    # Merge the logs so the sorting covers all log entries.
    merged_logs = itertools.chain(*logs)
    for entry in log_utils.filter_and_sort_log_entries(merged_logs):
        insert_id = entry['insertId']
        # skip duplicates
        if insert_id in seen_insert_ids:
            continue

        seen_insert_ids.add(insert_id)
        # Most of the pod name is the same across all pods, only the bit after the last '-' is unique.
        pod = entry['resource']['labels']['pod_name'].split('-')[-1]

        if log_utils.is_end_event(entry):
            if pod in in_progress_starts_by_pod:
                # Make a backpressure event object
                start = parser.isoparse(in_progress_starts_by_pod[pod])
                end = parser.isoparse(entry['timestamp'])
                event = BackpressureEvent(pod=pod, start=start, end=end)

                # Add this object to complete
                complete.append(event)

                # Remove the wip object from in_progress_starts_by_pod
                del in_progress_starts_by_pod[pod]

        elif log_utils.is_start_event(entry):
            # There are actually two timestamps in the JSON log entries which appear to represent different concepts:
            # time emitted ('jsonPayload.localTimestamp') versus time added to the log ('timestamp'). Time emitted would
            # seem to be preferable but that value is not specified with a timezone and is ambiguously interpreted by
            # the parsing code as being EST when it's actually UTC. This can make reading the report a bit confusing or
            # misleading. In practice the timestamps only seem to differ by small amounts, so no big deal to use
            # 'timestamp' with its explicit UTC timezone.
            in_progress_starts_by_pod[pod] = entry['timestamp']

    return complete
