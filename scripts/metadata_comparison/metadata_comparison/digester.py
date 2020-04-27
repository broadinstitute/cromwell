import argparse
import json
# import requests
from google.cloud import storage
import google.auth
# from googleapiclient.discovery import build as google_client_build
# import logging as log
import metadata_comparison.lib.argument_regex as reutil
import metadata_comparison.lib.util as util
# from datetime import timedelta, datetime, timezone
import dateutil.parser
import os
from pathlib import Path
from typing import Any, AnyStr, Dict, List, Union

Version = "0.0.1"
verbose = False


# Write the local and GCS code and then clean it up.

def main() -> None:
    args = parse_args()
    for path in args.local_paths:
        # GCS paths are currently returned as tuples of (bucket, object path) where local paths are just strings.
        # Useful, but a yucky implementation that should be cleaned up.
        if not isinstance(path, tuple):
            parent_path = Path(path)
            workflow_json_path = parent_path / 'workflow.json'
            # operations_dir = parent_path / 'operations'
            with open(workflow_json_path, 'r') as file:
                data = file.read()
                metadata = json.loads(data)
                digest_parent = parent_path / 'digests' / Version
                digest_path = digest_parent / 'digest.json'
                if not os.path.exists(digest_path) or args.force:
                    digest_parent.mkdir(parents=True, exist_ok=True)
                    with open(digest_path, 'w') as digest_file:
                        digested = digest(metadata)
                        json_string = json.dumps(digested, sort_keys=True, indent=4)
                        digest_file.write(json_string)
                else:
                    raise ValueError(f'digest file already exists at {digest_path} and --force not specified')
        else:
            gcs_bucket, gcs_object = path
            credentials, project_id = google.auth.default()
            storage_client = storage.Client(credentials=credentials)
            process_workflow_gcs(gcs_bucket, gcs_object, storage_client)


def process_workflow_gcs(gcs_bucket: AnyStr, gcs_path: AnyStr, storage_client: storage.Client) -> None:
    bucket = storage_client.get_bucket(gcs_bucket)
    blob = bucket.blob(f'{gcs_path}/extractor/metadata.json')
    json_string_bytes = blob.download_as_string()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Digest workflow metadata and job operation details, reading from and reuploading to GCS.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='whether to log verbosely (default False)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='whether to overwrite existing digests (default False)')
    parser.add_argument('local_paths', metavar="PATH", nargs='+', type=reutil.validate_gcs_if_gcs,
                        help="Location at which to find metadata (local or GCS)")

    return parser.parse_args()


JsonObject = Dict[str, Union[Union[None, AnyStr, float], Any]]
CallName = AnyStr


def digest(metadata: dict) -> dict:
    def call_fn(operation_mapping: Dict[CallName, JsonObject],
                operation_id: AnyStr,
                path: List[AnyStr],
                attempt: JsonObject) -> None:
        backend_status = attempt.get('backendStatus', 'Unknown')
        # This script should only ever be pointed at successful workflow metadata. All jobs that have a backend status
        # other than `Success` must have later been re-run successfully, so any un`Success`ful attempts are ignored.
        # It's possible that a future version of the digester might actually want to look at these jobs since they
        # may have completed some lifecycle events which could be useful in accumulating more performance data.
        if backend_status == 'Success':
            string_path = '.'.join(path)
            cromwell_start = attempt.get('start')
            cromwell_end = attempt.get('end')

            cromwell_total_time_seconds = (dateutil.parser.parse(cromwell_end) -
                                           dateutil.parser.parse(cromwell_start)).total_seconds()

            operation_mapping[string_path] = {
                "attempt": attempt.get('attempt'),
                "shardIndex": attempt.get('shardIndex'),
                "operationId": operation_id,
                "cromwellStart": cromwell_start,
                "cromwellEnd": cromwell_end,
                "cromwellTotalTimeSeconds": cromwell_total_time_seconds
            }

    shards = util.build_papi_operation_mapping(metadata, call_fn)
    return {'version': Version, 'calls': shards, 'workflowId': metadata['id']}


if __name__ == "__main__":
    main()
