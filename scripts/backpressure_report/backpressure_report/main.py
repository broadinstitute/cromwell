"""
Parse JSON logs for backpressure messages to produce a csv of backpressure seconds by time windows.
The Logs Explorer query looks like::

    resource.labels.container_name="cromwell1-runner-app"
    (jsonPayload.message=~"IoActor backpressure off" OR jsonPayload.message=~"Beginning IoActor backpressure")

Usage:

python main.py <files with json formatted Log Explorer logs>

"""
from datetime import timedelta
from dateutil import parser
import json
import sys

from lib.backpressure_event import BackpressureEvent


def build_log_jsons_from_input_files() -> list:
    return [json.load(open(f, 'r')) for f in sys.argv[1:]]


def is_event_start(entry) -> bool:
    return entry['jsonPayload']['message'].startswith('Beginning IoActor backpressure')


def is_event_end(entry) -> bool:
    return entry['jsonPayload']['message'] == 'IoActor backpressure off'


def filter_and_sort_log_entries(log) -> list:
    # Only start or end events are interesting
    filtered = [
        entry for entry in log if (is_event_start(entry) or is_event_end(entry))
    ]

    # Oldest first
    filtered.sort(key=(lambda e: e['timestamp']))
    return filtered


def build_backpressure_events_from_log_jsons(logs):
    complete = []
    seen_insert_ids = set(())
    in_progress_starts_by_pod = {}

    for log in logs:
        for entry in filter_and_sort_log_entries(log):
            insert_id = entry['insertId']
            # skip duplicates
            if insert_id in seen_insert_ids:
                continue

            seen_insert_ids.add(insert_id)
            pod = entry['resource']['labels']['pod_name'].split('-')[-1]

            if is_event_end(entry):
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

            elif is_event_start(entry):
                in_progress_starts_by_pod[pod] = entry['timestamp']

    return complete


def build_backpressure_windows_from_events(backpressure_events, window_width_in_hours=1):
    hour = backpressure_events[0].start.replace(minute=0, second=0, microsecond=0)
    next_hour = hour + timedelta(hours=window_width_in_hours)

    windows_by_hour = {hour: []}

    for window in backpressure_events:
        while window.start >= next_hour:
            hour = next_hour
            windows_by_hour[hour] = []
            next_hour = next_hour + timedelta(hours=window_width_in_hours)
        windows_by_hour[hour].append(window)
    return windows_by_hour


def print_windows(backpressure_windows):
    for interval, backpressure_events in backpressure_windows.items():
        print(f"{str(interval)},{sum([e.duration() for e in backpressure_events])}")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    log_jsons = build_log_jsons_from_input_files()
    events = build_backpressure_events_from_log_jsons(log_jsons)
    windows = build_backpressure_windows_from_events(events, window_width_in_hours=4)
    print_windows(windows)
