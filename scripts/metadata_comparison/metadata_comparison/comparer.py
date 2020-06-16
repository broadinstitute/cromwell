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
#   * pip3 install --upgrade pandas
#   * pip3 install --upgrade google-api-python-client
#   * pip3 install --upgrade google-cloud-storage
#
# Remember to login to create application default credentials before use:
#   % gcloud auth application-default login

from typing import AnyStr, List, Tuple
import argparse
import json
import pandas
import google.auth
from google.cloud import storage
import logging
from metadata_comparison.lib.comparison_paths import ComparisonPath, validate_path
from metadata_comparison.lib.logging import set_log_verbosity, quieten_chatty_imports
from metadata_comparison.lib.storage import upload_blob
from metadata_comparison.lib.argument_regex import gcs_path_regex_validator

logger = logging.getLogger('metadata_comparison.comparer')


def read_digester_jsons(path_strings: List[AnyStr],
                        _storage_client: storage.Client) -> List[Tuple[str, dict]]:

    result = []
    for path_string in path_strings:
        path = ComparisonPath.create(path_string)
        _json = json.loads(path.read_text())
        result.append((_json.get('workflowId'), _json))

    return result


def compare_jsons(_workflow_ids_and_jsons: List[Tuple[str, dict]]) -> pandas.DataFrame:
    """
    Uses pandas library to convert JSONs into dataframes, and concatenate those dataframes into a single one.
    Performs sanity check, producing exception, if at least one of the JSONs doesn't have matching subset of keys.
    """
    column_to_compare_name_ending = ".cromwellTotalTimeSeconds"
    version_column_name = "version"
    result = pandas.DataFrame()
    last_cols = []
    for workflow_id_and_json in _workflow_ids_and_jsons:
        df = pandas.json_normalize(workflow_id_and_json[1])
        cols = [c for c in df.columns if c.endswith(column_to_compare_name_ending)]
        cols.sort()
        cols.insert(0, version_column_name)

        if last_cols and last_cols != cols:
            at = _workflow_ids_and_jsons[0]
            msg = f"JSON data at {at} doesn't have matching subset of columns. Expected: {last_cols} but got {cols}"
            raise Exception(msg)

        last_cols = cols
        df.index = [workflow_id_and_json[0]]
        result = pandas.concat([result, df[cols]])

    rename_version_column_to = "digester format version"
    result.rename(columns={version_column_name: rename_version_column_to}, inplace=True)
    result.index.name = "workflow id"

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare performance metadata JSONs and produce CSV result')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('output_path', metavar='OUTPUT_PATH', type=validate_path, nargs=1,
                        help='Path to output CSV file.')
    parser.add_argument('digest_1', nargs=1, metavar='REQUIRED_DIGEST_1', type=validate_path,
                        help='First required digest path to compare.')
    parser.add_argument('digest_2', nargs=1, metavar='REQUIRED_DIGEST_2', type=validate_path,
                        help='Second required digest path to compare.')
    parser.add_argument('digests_more', nargs='*', metavar='OPTIONAL_MORE_DIGESTS', type=validate_path,
                        help='Any additional digest paths to compare.')

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()
    logger.info("Starting Comparer operation.")

    credentials, project_id = google.auth.default()
    storage_client = storage.Client(credentials=credentials)

    paths = args.digest_1 + args.digest_2 + args.digests_more
    workflow_ids_and_jsons = read_digester_jsons(paths, storage_client)
    comparison_result_df = compare_jsons(workflow_ids_and_jsons)
    result_csv_string = comparison_result_df.to_csv()

    out = ComparisonPath.create(args.output_path[0])
    out.write_text(result_csv_string)

    logger.info('Comparer operation completed successfully.')
