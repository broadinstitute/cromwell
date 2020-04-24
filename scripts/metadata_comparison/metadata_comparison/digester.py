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
import os
from pathlib import Path

Version = "0.0.1"
verbose = False


def main():
    args = parse_args()
    global verbose
    for path in args.local_paths:
        if not isinstance(path, tuple):
            parent_path = path if not path.endswith('/') else path[:-1]
            workflow_json_path = f'{parent_path}/workflow.json'
            # operations_dir = parent_path + 'operations'
            with open(workflow_json_path, 'r') as file:
                data = file.read()
                metadata = json.loads(data)
                digest_parent = Path(f'{parent_path}/digests/{Version}')
                digest_path = digest_parent / 'digest.json'
                if not os.path.exists(digest_path) or args.force:
                    digest_parent.mkdir(parents=True, exist_ok=True)
                    with open(digest_path, 'w') as digest_file:
                        digest_file.write(json.dumps(digest(metadata), sort_keys=True, indent=4))
                else:
                    raise ValueError(f'digest file already exists at {digest_path} and --force not specified')
        else:
            raise ValueError("Haven't written the GCS bits yet")


def parse_args():
    parser = argparse.ArgumentParser(
        description='Digest workflow metadata and job operation details, reading from and reuploading to GCS.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='whether to log verbosely (default False)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='whether to overwrite existing digests (default False)')
    parser.add_argument('local_paths', metavar="PATH", nargs='+', type=reutil.validate_gcs_if_gcs,
                        help="Location at which to find metadata (local or GCS)")

    return parser.parse_args()


def digest(metadata):
    def call_fn(operation_mapping, operation_id, path, attempt):
        backend_status = attempt.get('backendStatus', 'Unknown')
        if backend_status == 'Success':
            string_path = '.'.join(path)
            start = attempt.get('start')
            end = attempt.get('end')

            cromwell_total_time_seconds = (dateutil.parser.parse(end) - dateutil.parser.parse(start)).total_seconds()

            operation_mapping[string_path] = {
                "attempt": attempt.get('attempt'),
                "shardIndex": attempt.get('shardIndex'),
                "operationId": operation_id,
                "start": start,
                "end": end,
                "cromwellTotalTimeSeconds": cromwell_total_time_seconds
            }

    shards = util.build_papi_operation_mapping(metadata, call_fn)
    return { 'version': Version, 'calls': shards }


if __name__ == "__main__":
    main()
