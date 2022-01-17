import json


def build_log_jsons_from_input_files(input_files: list) -> list:
    # This could be written more compactly as a list comprehension, but I don't know how to close all the files and the
    # unit tests were rightly complaining about this.
    files = [open(f, 'r') for f in input_files]
    jsons = [json.load(f) for f in files]
    for f in files:
        f.close()
    return jsons


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
