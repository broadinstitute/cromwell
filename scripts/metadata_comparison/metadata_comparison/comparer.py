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

from typing import List, Tuple
import argparse
import json
import pandas
import google.auth
from google.cloud import storage
import logging
from metadata_comparison.lib.logging import set_log_verbosity, quieten_chatty_imports
from metadata_comparison.lib.storage import upload_blob
from metadata_comparison.lib.argument_regex import gcs_path_regex_validator, digester_version_regex_validator, \
    workflow_regex_validator

logger = logging.getLogger('metadata_comparison.comparer')

def read_digester_jsons_from_gcs(bucket_name: str,
                                 base_path: str,
                                 digester_version: str,
                                 workflow_ids: List[str],
                                 storage_client: storage.Client) -> List[Tuple[str, dict]]:
    bucket = storage_client.get_bucket(bucket_name)
    result = []
    for workflow_id in workflow_ids:
        blob = bucket.blob(f"{base_path}/{workflow_id}/digests/{digester_version}/digest.json")
        json_string_bytes = blob.download_as_string()
        result.append((workflow_id, json.loads(json_string_bytes)))

    return result


def compare_jsons(workflow_ids_and_jsons: List[Tuple[str, dict]]) -> pandas.DataFrame:
    """
    Uses pandas library to convert JSONs into dataframes, and concatenate those dataframes into a single one.
    Performs sanity check, producing exception, if at least one of the JSONs doesn't have matching subset of keys.
    """
    columnToCompareNameEnding = ".cromwellTotalTimeSeconds"
    versionColumnName = "version"
    result = pandas.DataFrame()
    last_cols = []
    for workflow_id_and_json in workflow_ids_and_jsons:
        df = pandas.json_normalize(workflow_id_and_json[1])
        cols = [c for c in df.columns if c.endswith(columnToCompareNameEnding)]
        cols.sort()
        cols.insert(0, versionColumnName)

        if last_cols and last_cols != cols:
            raise Exception(f"JSON data at {workflow_ids_and_jsons[0]} doesn't have matching subset of columns. Expected: {last_cols} but got {cols}")

        last_cols = cols
        df.index = [workflow_id_and_json[0]]
        result = pandas.concat([result, df[cols]])

    renameVersionColumnTo = "digester format version"
    result.rename(columns={versionColumnName: renameVersionColumnTo}, inplace=True)
    result.index.name = "workflow id"

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare performance metadata JSONs and produce CSV result')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--digester-version', metavar='DIGESTERVERSION', type=digester_version_regex_validator, nargs=1,
                        help='Compare digests produced by this version of the digester')
    parser.add_argument('--digest-gcs-base-path', metavar='DIGESTGCSBASEPATH', type=gcs_path_regex_validator, nargs=1,
                        help='GCS base path to the directory containing JSONs produced by digester')
    parser.add_argument('--output-gcs-file-path', metavar='OUTPUTGCSFILE', type=gcs_path_regex_validator, nargs=1,
                        help='GCS path to output CSV file')
    parser.add_argument('--workflow-ids', metavar='WORKFLOWIDS', type=workflow_regex_validator, nargs='+',
                        help='Workflow ids for performance comparison')

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()
    logger.info("Starting Comparer operation.")

    credentials, project_id = google.auth.default()
    storage_client = storage.Client(credentials = credentials)
    input_gcs_bucket, input_gcs_path = args.digest_gcs_base_path[0]

    workflow_ids_and_jsons = read_digester_jsons_from_gcs(input_gcs_bucket, input_gcs_path, args.digester_version[0], args.workflow_ids, storage_client)
    comparison_result_df = compare_jsons(workflow_ids_and_jsons)
    result_csv_string = comparison_result_df.to_csv()

    output_gcs_bucket, output_gcs_path = args.output_gcs_file_path[0]
    upload_blob(output_gcs_bucket, result_csv_string, output_gcs_path, storage_client, logger)

    logger.info('Comparer operation completed successfully.')
