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


def read_input_files_to_log_jsons():
    return [json.load(open(f, 'r')) for f in sys.argv[1:]]


def parse_log_jsons_to_backpressure_events(logs):
    """
    Logs are assumed to be sorted *most recent first*, this will not work as is if sorted any other way.
    :param logs:
    :return:
    """
    complete = []
    seen_insert_ids = set(())
    in_progress_ends_by_pod = {}

    for log in logs:
        for entry in log:
            insert_id = entry['insertId']
            # skip duplicates
            if insert_id in seen_insert_ids:
                continue
            else:
                seen_insert_ids.add(insert_id)

            pod = entry['resource']['labels']['pod_name'].split('-')[-1]
            message = entry['jsonPayload']['message']

            if message.startswith('Beginning IoActor backpressure'):
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

            elif message.startswith('IoActor backpressure off'):
                in_progress_ends_by_pod[pod] = entry['timestamp']

    return complete


def build_windows_by_hour(windows, window_width_hours=1):
    reversed_windows = windows.copy()
    reversed_windows.reverse()
    hour = reversed_windows[0].start.replace(minute=0, second=0, microsecond=0)
    next_hour = hour + timedelta(hours=window_width_hours)

    windows_by_hour = {hour: []}

    for window in reversed_windows:
        while window.start >= next_hour:
            hour = next_hour
            windows_by_hour[hour] = []
            next_hour = next_hour + timedelta(hours=window_width_hours)
        windows_by_hour[hour].append(window)
    return windows_by_hour


def print_windows(by_hour):
    for hour, events in by_hour.items():
        print(str(hour), sum([e.duration() for e in events]))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    json_logs = read_input_files_to_log_jsons()
    events = parse_log_jsons_to_backpressure_events(json_logs)
    by_hour = build_windows_by_hour(events, window_width_hours=4)
    print_windows(by_hour)
