from backpressure_report.lib import log_utils
from datetime import datetime
from dateutil import parser
from typing import AnyStr


class BackpressureEvent:
    def __init__(self, pod: str, start: datetime, end: datetime):
        self.pod = pod
        self.start = start
        self.end = end

    def duration(self):
        return (self.end - self.start).seconds

    def __str__(self) -> AnyStr:
        return f"BackpressureEvent(pod = {self.pod},start={str(self.start)},duration={(self.duration())}s)"


def build_backpressure_events_from_log_jsons(logs):
    complete = []
    seen_insert_ids = set(())
    in_progress_starts_by_pod = {}

    for log in logs:
        for entry in log_utils.filter_and_sort_log_entries(log):
            insert_id = entry['insertId']
            # skip duplicates
            if insert_id in seen_insert_ids:
                continue

            seen_insert_ids.add(insert_id)
            pod = entry['resource']['labels']['pod_name'].split('-')[-1]

            if log_utils.is_event_end(entry):
                if pod in in_progress_starts_by_pod.keys():
                    # Make a backpressure event object
                    start = parser.isoparse(in_progress_starts_by_pod[pod])
                    end = parser.isoparse(entry['timestamp'])
                    event = BackpressureEvent(pod=pod, start=start, end=end)

                    # Add this object to complete
                    complete.append(event)

                    # Remove the wip object from in_progress_starts_by_pod
                    in_progress_starts_by_pod.pop(pod)
                    # print(event)

            elif log_utils.is_event_start(entry):
                in_progress_starts_by_pod[pod] = entry['timestamp']

    return complete
