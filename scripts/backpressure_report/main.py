"""
Parse a time-sorted array of JSON log entries for backpressure messages. The Logs Explorer query looks like:

resource.labels.container_name="cromwell1-runner-app"
(jsonPayload.message=~"IoActor backpressure off" OR jsonPayload.message=~"Beginning IoActor backpressure")

"""
from datetime import timedelta
from dateutil import parser
import json
import sys

from lib.backpressure_event import BackpressureEvent


def build_log_jsons_from_input_files():
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

    # Most recent first
    filtered.sort(key=(lambda e: e['timestamp']), reverse=True)
    return filtered


def build_backpressure_events_from_log_jsons(logs):
    complete = []
    seen_insert_ids = set(())
    in_progress_ends_by_pod = {}

    for log in logs:
        for entry in filter_and_sort_log_entries(log):
            insert_id = entry['insertId']
            # skip duplicates
            if insert_id in seen_insert_ids:
                continue

            seen_insert_ids.add(insert_id)
            pod = entry['resource']['labels']['pod_name'].split('-')[-1]

            if is_event_start(entry):
                if pod in in_progress_ends_by_pod.keys():
                    # Make a backpressure event object
                    end = parser.isoparse(in_progress_ends_by_pod[pod])
                    start = parser.isoparse(entry['timestamp'])
                    event = BackpressureEvent(pod=pod, start=start, end=end)

                    # Add this object to complete
                    complete.append(event)

                    # Remove the wip object from in_progress_ends_by_pod
                    in_progress_ends_by_pod.pop(pod)
                    # print(event)

            elif is_event_end(entry):
                in_progress_ends_by_pod[pod] = entry['timestamp']

    return complete


def build_backpressure_windows_from_events(windows, window_width_in_hours=1):
    reversed_windows = windows.copy()
    reversed_windows.reverse()
    hour = reversed_windows[0].start.replace(minute=0, second=0, microsecond=0)
    next_hour = hour + timedelta(hours=window_width_in_hours)

    windows_by_hour = {hour: []}

    for window in reversed_windows:
        while window.start >= next_hour:
            hour = next_hour
            windows_by_hour[hour] = []
            next_hour = next_hour + timedelta(hours=window_width_in_hours)
        windows_by_hour[hour].append(window)
    return windows_by_hour


def print_windows(by_hour):
    for hour, events in by_hour.items():
        print(str(hour), sum([e.duration() for e in events]))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    log_jsons = build_log_jsons_from_input_files()
    events = build_backpressure_events_from_log_jsons(log_jsons)
    by_hour = build_backpressure_windows_from_events(events, window_width_in_hours=4)
    print_windows(by_hour)
