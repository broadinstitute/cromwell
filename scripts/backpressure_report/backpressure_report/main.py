from backpressure_report.lib import backpressure_event, backpressure_window, log_utils
import sys


if __name__ == '__main__':
    log_jsons = log_utils.build_log_jsons_from_input_files(sys.argv[1:])
    events = backpressure_event.build_backpressure_events_from_log_jsons(log_jsons)
    windows = backpressure_window.build_windows_from_events(events, window_width_in_hours=1)
    backpressure_window.print_windows(windows)