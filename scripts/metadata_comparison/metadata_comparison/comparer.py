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

logger = logging.getLogger('metadata_comparison.comparer')


def read_digester_jsons(path_strings: List[AnyStr],
                        _storage_client: storage.Client) -> List[Tuple[str, dict]]:

    result = []
    for path_string in path_strings:
        path = ComparisonPath.create(path_string)
        _json = json.loads(path.read_text())
        result.append((_json.get('workflowId'), _json))

    return result


def compare_jsons() -> List[List[AnyStr]]:
    columns = [
        {'PAPI Total Time': '.papiTotalTimeSeconds'}
    ]
    return [["foo"]]


def check_digest_versions():
    version_1, version_2 = [j.get('version') for j in [json_1, json_2]]

    if version_1 != version_2:
        msg = f"Inconsistent digest versions: First JSON digest is {version_1} but second is {version_2}"
        raise ValueError(msg)

    call_keys_1, call_keys_2 = [j.get('calls').keys() for j in [json_1, json_2]]

    if call_keys_1 != call_keys_2:
        raise ValueError("Digest files do not have the same call keys, not the same workflow and/or sample.")


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
    check_digest_versions()

    comparison_data = compare_jsons()
    out = ComparisonPath.create(args.output_path[0])
    out.write_text(csv_string_from_data(comparison_data))

    logger.info('Comparer operation completed successfully.')
