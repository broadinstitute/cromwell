import argparse
import json
import requests
from google.cloud import storage
import google.auth
from googleapiclient.discovery import build as google_client_build
import logging as log
import metadata_comparison.lib.argument_regex as reutil
import metadata_comparison.lib.util as util
from datetime import timedelta, datetime, timezone
import dateutil.parser

Version = "0.0.1"
verbose = False


def main():
    args = parse_args()
    global verbose
    verbose = args.verbose
    print(args)
    for gcs_path in args.gcs_paths:
        parent_path = f'gs://{gcs_path[0]}/{gcs_path[1]}/'
        workflow_json_path = parent_path + 'workflow.json'
        operations_path = parent_path + 'operations'
        print(f'would expect to find a workflow json at {workflow_json_path} and operations under {operations_path}')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Digest workflow metadata and job operation details, reading from and reuploading to GCS.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='whether to log verbosely (default False)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='whether to overwrite existing digests (default False)')
    parser.add_argument('gcs_paths', metavar='GCSPATH', type=reutil.gcs_path_regex_validator, nargs='+',
                        help='GCS paths containing workflow data')

    return parser.parse_args()


def find_succeeded_shards(metadata):
    def call_fn(operation_mapping, operation_id, path, attempt):
        backend_status = attempt.get('backendStatus', 'Unknown')
        # print(f"digester seeing path of {path} and backendStatus {backend_status}")
        if backend_status == 'Success':
            # print("Job was Succeeded woot")
            string_path = '.'.join(path)
            start = attempt.get('start')
            end = attempt.get('end')

            overall_runtime_seconds = (dateutil.parser.parse(end) - dateutil.parser.parse(start)).total_seconds()

            operation_mapping[string_path] = {
                "attempt": attempt.get('attempt'),
                "shardIndex": attempt.get('shardIndex'),
                "operationId": operation_id,
                "start": start,
                "end": end,
                "overallRuntimeSeconds": overall_runtime_seconds
            }

    shards = util.build_papi_operation_mapping(metadata, call_fn)
    print(shards)


if __name__ == "__main__":
    main()
