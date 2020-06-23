#!/usr/bin/env python3
#
# comparer.py
#
# Purpose: Compare performance metadata JSON files produced by Digester and produce result in CSV format
#
# Usage: python3 comparer.py [-h] [-v] [--json_paths JSONPATH [JSONPATH ...]] [--output_path OUTPUTPATH]
#
# Python Prereqs (at least, the ones which I needed to manually install... YMMV):
#
#   * pip3 install --upgrade google-api-python-client
#   * pip3 install --upgrade google-cloud-storage
#
# Remember to login to create application default credentials before use:
#   % gcloud auth application-default login

import argparse
import json
import logging
from metadata_comparison.lib.digester_keys import *
from metadata_comparison.lib.comparison_paths import ComparisonPath, validate_path
from metadata_comparison.lib.logging import set_log_verbosity, quieten_chatty_imports
from metadata_comparison.lib.operation_ids import JsonObject
from typing import AnyStr, List

logger = logging.getLogger('metadata_comparison.comparer')


class DigesterKey:
    def __init__(self, json_key: AnyStr, display_text: AnyStr):
        self.json_key = json_key
        self.display_text = display_text


DigesterKeys = [
    DigesterKey(json_key=PapiTotalTimeSeconds, display_text='Total PAPI time (seconds)'),
    DigesterKey(json_key=StartupTimeSeconds, display_text='Startup (seconds)'),
    DigesterKey(json_key=DockerImagePullTimeSeconds, display_text='Docker Pull (seconds)'),
    DigesterKey(json_key=LocalizationTimeSeconds, display_text='Localization (seconds)'),
    DigesterKey(json_key=UserCommandTimeSeconds, display_text='User command (seconds)'),
    DigesterKey(json_key=DelocalizationTimeSeconds, display_text='Delocalization (seconds)'),
    DigesterKey(json_key=OtherTimeSeconds, display_text='Other time (seconds)'),
    DigesterKey(json_key=MachineType, display_text='Machine type')
]


MachineTypes = {
    'n1-highcpu-16': 	0.5672,
    'n1-highmem-2': 	0.12,
    'n1-standard-1':	0.05,
    'n1-standard-2':	0.10,
    'g1-small':	0.03
}


class CallKey:
    def __init__(self, full, call_prefixes: List[AnyStr]):
        self.full = full
        # Call keys without repetitive prefixes are useful for building comparer reports.
        self.without_prefix = full
        for prefix in call_prefixes:
            if self.full.startswith(prefix):
                self.without_prefix = self.full[len(prefix):]
                break


def format_seconds(total_seconds: float) -> AnyStr:
    hours, remainder = divmod(int(total_seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours}:{minutes:02}:{seconds:02}'


def compare_jsons(json_1: JsonObject, json_2: JsonObject,
                  name_1: AnyStr, name_2: AnyStr,
                  call_prefixes_to_remove: List[AnyStr]) -> List[List[AnyStr]]:
    call_keys = [CallKey(k, call_prefixes_to_remove) for k in json_1.get('calls').keys()]

    def without_prefix(call_key: CallKey) -> AnyStr:
        return call_key.without_prefix

    call_keys_sorted_without_prefix = sorted(call_keys, key=without_prefix)

    rows = []
    header_row = []
    percent_contrib_row = []
    total_row = []
    machine_type_weighted_row = []
    first = True
    for call_key in call_keys_sorted_without_prefix:
        row = []
        call_1, call_2 = [j.get('calls').get(call_key.full) for j in [json_1, json_2]]
        if first:
            header_row.append('job')
            percent_contrib_row.append('% contribution to total run time')
            total_row.append('Total')
            machine_type_weighted_row.append('Machine-type weighted total (cost-ish)')
        row.append(call_key.without_prefix)

        if first:
            header_row.append('Machine type')
            percent_contrib_row.append('')
            total_row.append('')
            machine_type_weighted_row.append('')
        row.append(call_1.get(MachineType))

        time_1 = call_1.get(PapiTotalTimeSeconds)
        time_2 = call_2.get(PapiTotalTimeSeconds)
        row.append(format_seconds(time_1))
        row.append(format_seconds(time_2))
        row.append(f'{((time_2 - time_1) / time_1) * 100:.2f}')
        if first:
            for name in [name_1, name_2]:
                header_row.append(f'{name} Total PAPI Time')
            header_row.append(f'Percent increase Total PAPI time {name_1} -> {name_2}')
            percent_contrib_row.append(f'{100:.2f}')

            def sum_call_times(j: JsonObject, key: AnyStr):
                times = [j.get('calls').get(ck.full).get(key) for ck in call_keys]
                return sum(times)

            total_time_1, total_time_2 = [
                sum_call_times(j, PapiTotalTimeSeconds) for j in [json_1, json_2]]

            total_percent_2 = (total_time_2 / total_time_1) * 100
            percent_contrib_row.append(f'{total_percent_2:.2f}')
            total_percent_percent = ((total_time_2 - total_time_1) * 100) / total_time_1
            percent_contrib_row.append(f'{total_percent_percent:.2f}')

            total_row.append(format_seconds(total_time_1))
            total_row.append(format_seconds(total_time_2))
            total_row.append(f'{total_percent_percent:.2f}')

        rows.append(row)

        first = False

    rows.insert(0, machine_type_weighted_row)
    rows.insert(0, total_row)
    rows.insert(0, percent_contrib_row)
    rows.insert(0, header_row)
    return rows


def error_checks(json_1: JsonObject, json_2: JsonObject):
    version_1, version_2 = [j.get('version') for j in [json_1, json_2]]

    if version_1 != version_2:
        msg = f"Inconsistent digest versions: First JSON digest is {version_1} but second is {version_2}"
        raise ValueError(msg)

    call_keys_1, call_keys_2 = [list(j.get('calls').keys()) for j in [json_1, json_2]]
    call_keys_1.sort()
    call_keys_2.sort()

    if call_keys_1 != call_keys_2:
        raise ValueError('The specified digest files do not have the same call keys. These digests cannot be ' +
                         'compared and probably are not from the same workflow and sample.')

    def machine_types(j: JsonObject, keys: List[AnyStr]) -> List[AnyStr]:
        return [j.get('calls').get(k).get('machineType') for k in keys]

    machine_types_1, machine_types_2 = [machine_types(j, call_keys_1) for j in [json_1, json_2]]

    if machine_types_1 != machine_types_2:
        # This should perhaps become a warning or something the user can --force their way through if we want to
        # allow experimentation with different machine types to reduce cost, but for the time being it's unexpected.
        raise ValueError('The specified digest files cannot be meaningfully compared as they contain calls with '
                         'different machine types for corresponding jobs.')

    for call_key in call_keys_1:
        call_1 = json_1.get('calls').get(call_key)
        call_2 = json_2.get('calls').get(call_key)
        for call in [call_1, call_2]:
            for digester_key in DigesterKeys:
                json_key = digester_key.json_key
                if json_key not in call:
                    if call == call_1:
                        nth = "first"
                        json_file = args.digest1[0]
                    else:
                        nth = "second"
                        json_file = args.digest2[0]
                    raise ValueError(
                        f"In {nth} digest JSON '{json_file}': call '{call_key}' does not contain required key '{json_key}'")


def json_from_path_string(path_string: AnyStr) -> JsonObject:
    path = ComparisonPath.create(path_string)
    return json.loads(path.read_text())


def csv_string_from_data(data: List[List[AnyStr]]) -> AnyStr:
    rows = [','.join(row) for row in data]
    return '\n'.join(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare performance digest JSONs and produce CSV result')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--name1', nargs=1, metavar='NAME_FOR_DIGEST_1', type=str,
                        required=True, help='Name to use for the first digest')
    parser.add_argument('--name2', nargs=1, metavar='NAME_FOR_DIGEST_2', type=str,
                        required=True, help='Name to use for the second digest')
    parser.add_argument('--digest1', nargs=1, metavar='PATH_TO_DIGEST_1', type=validate_path,
                        required=True, help='First digest path to compare.')
    parser.add_argument('--digest2', nargs=1, metavar='PATH_TO_DIGEST_2', type=validate_path,
                        required=True, help='Second digest path to compare.')
    parser.add_argument('--output-path', metavar='OUTPUT_PATH', type=validate_path, nargs=1,
                        required=True, help='Path to output CSV file.')
    parser.add_argument('--call-prefix-to-remove', metavar='CALL_PREFIX_TO_REMOVE', type=str, nargs='*',
                        help='Call prefix to remove if present.')

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()
    logger.info("Starting Comparer operation.")

    _json_1, _json_2 = [json_from_path_string(p[0]) for p in [args.digest1, args.digest2]]
    error_checks(_json_1, _json_2)

    prefixes = [] if not args.call_prefix_to_remove else args.call_prefix_to_remove

    comparison_data = compare_jsons(_json_1, _json_2, args.name1[0], args.name2[0], prefixes)
    out = ComparisonPath.create(args.output_path[0])
    out.write_text(csv_string_from_data(comparison_data))

    logger.info('Comparer operation completed successfully.')
