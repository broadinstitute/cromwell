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

from typing import AnyStr, List, Tuple
import argparse
import json
from google.cloud import storage
import logging
from metadata_comparison.lib.comparison_paths import ComparisonPath, validate_path
from metadata_comparison.lib.logging import set_log_verbosity, quieten_chatty_imports
from metadata_comparison.lib.operation_ids import JsonObject
from metadata_comparison.lib.digester_keys import *

logger = logging.getLogger('metadata_comparison.comparer')

DigesterKeys = [
    {
        'jsonKey': PapiTotalTimeSeconds,
        'displayText': "Total PAPI time (seconds)"
    },
    {
        'jsonKey': StartupTimeSeconds,
        'displayText': 'Startup (seconds)'
    },
    {
        'jsonKey': DockerImagePullTimeSeconds,
        'displayText': 'Docker Pull (seconds)'
    },
    {
        'jsonKey': LocalizationTimeSeconds,
        'displayText': 'Localization (seconds)'
    },
    {
        'jsonKey': UserCommandTimeSeconds,
        'displayText': 'User command (seconds)'
    },
    {
        'jsonKey': DelocalizationTimeSeconds,
        'displayText': 'Delocalization (seconds)'
    },
    {
        'jsonKey': OtherTimeSeconds,
        'displayText': 'Other time (seconds)'
    },
]


def read_digester_jsons(path_strings: List[AnyStr],
                        _storage_client: storage.Client) -> List[Tuple[str, dict]]:
    result = []
    for path_string in path_strings:
        path = ComparisonPath.create(path_string)
        _json = json.loads(path.read_text())
        result.append((_json.get('workflowId'), _json))

    return result


def compare_jsons() -> List[List[AnyStr]]:
    return [["foo"]]


def error_checks():
    version_1, version_2 = [j.get('version') for j in [json_1, json_2]]

    if version_1 != version_2:
        msg = f"Inconsistent digest versions: First JSON digest is {version_1} but second is {version_2}"
        raise ValueError(msg)

    call_keys_1, call_keys_2 = [list(j.get('calls').keys()) for j in [json_1, json_2]]
    call_keys_1.sort()
    call_keys_2.sort()

    if call_keys_1 != call_keys_2:
        raise ValueError('The specified digest files do not have the same call keys. These digests cannot be ' +
                         'compared and likely do not derive from the same workflow and/or sample.')

    def machine_types(j: JsonObject, keys: List[AnyStr]) -> List[AnyStr]:
        return [j.get('calls').get(k).get('machineType') for k in keys]

    machine_types_1, machine_types_2 = [machine_types(j, call_keys_1) for j in [json_1, json_2]]

    if machine_types_1 != machine_types_2:
        raise ValueError('The specified digest files cannot be meaningfully compared as they contain calls with '
                         'different machine types for corresponding jobs.')

    for call_key in call_keys_1:
        call_1 = json_1.get('calls').get(call_key)
        call_2 = json_2.get('calls').get(call_key)
        for call in [call_1, call_2]:
            for digester_key in DigesterKeys:
                json_key = digester_key.get('jsonKey')
                if json_key not in call_1:
                    if call == call_1:
                        nth = "first"
                        json_file = args.digest1
                    else:
                        nth = "second"
                        json_file = args.digest2
                    raise ValueError(
                        f"In {nth} JSON '{json_file[0]}': call '{call_key}' does not contain required key '{json_key}'")


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

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()
    logger.info("Starting Comparer operation.")

    json_1, json_2 = [json_from_path_string(p[0]) for p in [args.digest1, args.digest2]]
    name_1, name_2 = args.name1, args.name2
    error_checks()

    comparison_data = compare_jsons()
    out = ComparisonPath.create(args.output_path[0])
    out.write_text(csv_string_from_data(comparison_data))

    logger.info('Comparer operation completed successfully.')
