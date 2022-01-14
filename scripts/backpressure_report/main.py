"""
Parse a time-sorted array of JSON log entries for backpressure messages. The Logs Explorer query looks like:

resource.labels.container_name="cromwell1-runner-app"
(jsonPayload.message=~"IoActor backpressure off" OR jsonPayload.message=~"Beginning IoActor backpressure")

And the jq cleanup code looks like:

jq '[
  .[] |
  {
    message: .jsonPayload.message,
    timestamp: .jsonPayload.localTimestamp,
    pod:.resource.labels.pod_name
  }] |
reverse' downloaded-logs-20220112-191349.json > logs-20220112-191349-processed.json

"""
import argparse
from collections import defaultdict
from datetime import timedelta
from dateutil import parser
import json


def parse_args():
    """
    Parse the args, for now only defines 'input' as the name of the input file.
    :return: args dictionary
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("input")
    args = parser.parse_args()
    return args


def read_input(infile):
    """
    Read the specified input file as return as JSON.
    :param infile:
    :return:
    """
    f = open(infile, 'r')
    return json.load(f)


def parse_logs_to_backpressure_windows(logs):
    """
    Parse the processed logs to a list of backpressure window sorted by start:
    {
      "start": <start timestamp>,
      "end": <end timestamp>,
      "pod": <pod name>
    }
    :param logs:
    :return:
    """

    in_progress = {}
    complete = []

    for entry in logs:
        pod = entry['pod'].split('-')[-1]
        message = entry['message']
        if message == 'IoActor backpressure off':
            if pod in in_progress:
                # Make a backpressure window object
                start = parser.isoparse(in_progress[pod])
                end = parser.isoparse(entry['timestamp'])

                obj = {
                    'start': start,
                    'end': end,
                    'pod': pod,
                    'duration': (end - start).seconds
                }
                # Add this object to complete
                complete.append(obj)
                # Remove the wip object from in_progress
                in_progress.pop(pod)
                # print(obj)
        elif message.startswith('Beginning IoActor backpressure'):
            in_progress[pod] = entry['timestamp']
    return complete


def build_windows_by_hour(windows, hour_delta=1):
    hour = windows[0]['start'].replace(minute=0, second=0, microsecond=0)
    next_hour = hour + timedelta(hours=hour_delta)

    windows_by_hour = {hour: []}

    for window in windows:
        while window['start'] >= next_hour:
            hour = next_hour
            windows_by_hour[hour] = []
            next_hour = next_hour + timedelta(hours=hour_delta)
        windows_by_hour[hour].append(window)
    return windows_by_hour


def print_windows_by_hour(by_hour):
    for hour, windows in by_hour.items():
        print(str(hour), sum([w['duration'] for w in windows]))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    args = parse_args()
    json_logs = read_input(args.input)
    windows = parse_logs_to_backpressure_windows(json_logs)
    by_hour = build_windows_by_hour(windows, hour_delta=24)
    print_windows_by_hour(by_hour)
