"""
Parse JSON logs for backpressure messages to produce a csv of backpressure seconds by time windows.
The Logs Explorer query looks like::

    resource.labels.container_name="cromwell1-runner-app"
    (jsonPayload.message=~"IoActor backpressure off" OR jsonPayload.message=~"Beginning IoActor backpressure")

Usage:

python -m backpressure_report.main. <files with json formatted Log Explorer logs>

"""
from backpressure_report.lib import backpressure_event, log_utils
from datetime import timedelta
import sys


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
    log_jsons = log_utils.build_log_jsons_from_input_files(sys.argv[1:])
    events = backpressure_event.build_backpressure_events_from_log_jsons(log_jsons)
    windows = build_backpressure_windows_from_events(events, window_width_in_hours=4)
    print_windows(windows)
