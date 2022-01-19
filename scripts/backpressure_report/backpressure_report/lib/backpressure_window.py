from collections import defaultdict
from datetime import timedelta


class BackpressureWindow:
    def __init__(self, timestamp):
        self.timestamp = timestamp
        self.pod_events = defaultdict(list)

    def add_event(self, event):
        self.pod_events[event.pod].append(event)

    def durations_by_pod(self) -> dict:
        # defaultdict(int) also returns 0 by default, more cryptically
        ret = defaultdict(lambda _: 0)
        for pod, events in self.pod_events.items():
            ret[pod] = sum([e.duration() for e in events])
        return ret


def build_windows_from_events(backpressure_events, window_width_in_hours=1):
    """
    Generate barchart-friendly time windows with counts of backpressuring durations within each window.

    :param backpressure_events: a list of BackpressureEvents to be broken up into time windows
    :param window_width_in_hours: how wide each time window should be in hours
    :return: a dictionary with timestamp keys to list of BackpressureEvent values
    """

    # The logic below is highly dependent on events being sorted by start timestamp oldest to newest.
    sorted_events = backpressure_events.copy()
    sorted_events.sort(key=lambda e: e.start)

    interval = sorted_events[0].start.replace(minute=0, second=0, microsecond=0)
    next_interval = interval + timedelta(hours=window_width_in_hours)

    all_pods = set(())
    windows = [BackpressureWindow(interval)]

    for event in sorted_events:
        all_pods.add(event.pod)
        while event.start >= next_interval:
            interval = next_interval
            windows.append(BackpressureWindow(interval))
            next_interval = next_interval + timedelta(hours=window_width_in_hours)
        windows[-1].add_event(event)
    return windows, all_pods


def print_windows(windows):
    """
    CSV format output generation for the specified backpressure windows.
    :param windows: dictionary of time windows to list of BackpressureEvents
    """
    for interval, backpressure_events in windows.items():
        print(f"{str(interval)},{sum([e.duration() for e in backpressure_events])}")
