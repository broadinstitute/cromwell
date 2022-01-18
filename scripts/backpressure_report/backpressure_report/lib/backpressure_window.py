from datetime import timedelta


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

    hour = sorted_events[0].start.replace(minute=0, second=0, microsecond=0)
    next_hour = hour + timedelta(hours=window_width_in_hours)

    windows_by_hour = {hour: []}

    for window in sorted_events:
        while window.start >= next_hour:
            hour = next_hour
            windows_by_hour[hour] = []
            next_hour = next_hour + timedelta(hours=window_width_in_hours)
        windows_by_hour[hour].append(window)
    return windows_by_hour


def print_windows(windows):
    """
    CSV format output generation for the specified backpressure windows.
    :param windows: dictionary of time windows to list of BackpressureEvents
    """
    for interval, backpressure_events in windows.items():
        print(f"{str(interval)},{sum([e.duration() for e in backpressure_events])}")
