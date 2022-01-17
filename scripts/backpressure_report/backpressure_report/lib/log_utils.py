import json


def build_log_jsons_from_input_files(input_files: list) -> list:
    def load(path):
        # `with` to auto-close files
        with open(path, 'r') as f:
            return json.load(f)

    return [load(f) for f in input_files]


def is_start_event(entry) -> bool:
    return entry['jsonPayload']['message'].startswith('Beginning IoActor backpressure')


def is_end_event(entry) -> bool:
    return entry['jsonPayload']['message'] == 'IoActor backpressure off'


def filter_and_sort_log_entries(log) -> list:
    # Only start or end events are interesting
    filtered = [
        entry for entry in log if (is_start_event(entry) or is_end_event(entry))
    ]

    # Oldest first
    filtered.sort(key=(lambda e: e['timestamp']))
    return filtered
