#!/usr/bin/env python3
#
# extractor.py
#
# Purpose: Read workflow metadata from Cromwell, and all metadata for its jobs,
# and upload it to a GCS bucket
#
# Usage: python3 extractor.py <GCS path> <workflowId> [<workflowId2> [...]]
#
# Python Prereqs (at least, the ones which I needed to manually install... YMMV):
#
#   * pip3 install --upgrade requests
#   * pip3 install --upgrade google-api-python-client
#   * pip3 install --upgrade google-cloud
#   * pip3 install --upgrade google-cloud-storage
#
# Remember to login to create application default credentials before use:
#   % gcloud auth application-default login

import argparse
import json
import requests
from google.cloud import storage
import google.auth
import logging
from metadata_comparison.lib.argument_regex import url_regex_validator, gcs_path_regex_validator, workflow_regex_validator
from metadata_comparison.lib.operation_ids import get_operation_id_number, find_operation_ids_in_metadata
from metadata_comparison.lib.papi.papi_clients import PapiClients

logger = logging.getLogger('metadata_comparison.extractor')


def set_log_verbosity(verbose: bool) -> None:
    if verbose:
        logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=logging.WARNING)


def quieten_chatty_imports() -> None:
    logging.getLogger('googleapiclient.discovery_cache').setLevel(logging.ERROR)
    logging.getLogger('googleapiclient.discovery').setLevel(logging.WARNING)


def upload_local_checkout():
    #   # Make a snapshot of the local Cromwell git repo
    #
    #   hash=$(git rev-parse HEAD)
    #
    #   # Make sure all files are checked in.
    #   if [[ $(git status --short | wc -l) -ne 0 ]]; then echo "stuff needs to be checked in"; exit 1; fi
    #
    #   git_snapshots_dir="../git-snapshots"
    #   mkdir -p "${git_snapshots_dir}"
    #   git_snapshot="${git_snapshots_dir}/$hash"
    #   git clone . --local "${git_snapshot}"
    raise Exception("Not Implemented")


def fetch_raw_workflow_metadata(cromwell_url: str, workflow: str) -> (requests.Response, dict):
    """Fetches workflow metadata for a workflow. Returns the raw response and the dict read from json"""
    url = f'{cromwell_url}/api/workflows/v1/{workflow}/metadata?expandSubWorkflows=true'
    logger.info(f'Fetching Cromwell metadata from {url}...')
    result = requests.get(url)
    return result.content, result.json()


def upload_blob(bucket_name: str, source_file_contents: bytes, destination_blob_name: str, gcs_storage_client: storage.Client) -> None:
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"

    bucket = gcs_storage_client.bucket(bucket_name)

    logger.info(f'Uploading file content to gs://{bucket_name}/{destination_blob_name}...')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)


def upload_workflow_metadata_json(bucket_name: str, raw_workflow_metadata: bytes, workflow_gcs_base_path: str, gcs_storage_client: storage.Client) -> None:
    workflow_gcs_metadata_upload_path = f'{workflow_gcs_base_path}/metadata.json'
    upload_blob(bucket_name, raw_workflow_metadata, workflow_gcs_metadata_upload_path, gcs_storage_client)


def upload_operations_metadata_json(bucket_name: str, operation_id: str, operations_metadata: dict, workflow_gcs_base_path: str, gcs_storage_client: storage.Client) -> None:
    """Uploads metadata to cloud storage, as json"""
    operation_upload_path = f'{workflow_gcs_base_path}/operations/{get_operation_id_number(operation_id)}.json'
    formatted_metadata = json.dumps(operations_metadata, indent=2)
    upload_blob(bucket_name, bytes(formatted_metadata, 'utf-8'), operation_upload_path, gcs_storage_client)


def process_workflow(cromwell_url: str, gcs_bucket: str, gcs_path: str, gcs_storage_client: storage.Client, papi_clients: PapiClients, workflow: str) -> None:
    raw_metadata, json_metadata = fetch_raw_workflow_metadata(cromwell_url, workflow)
    workflow_gcs_base_path = f'{gcs_path}/{workflow}/extractor'

    operation_ids = find_operation_ids_in_metadata(json_metadata)
    for id in operation_ids:
        operation_metadata = papi_clients.request_operation_metadata(id)
        upload_operations_metadata_json(gcs_bucket, id, operation_metadata, workflow_gcs_base_path, gcs_storage_client)
    upload_workflow_metadata_json(gcs_bucket, raw_metadata, workflow_gcs_base_path, gcs_storage_client)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract metadata and operation details for workflows and upload to GCS')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('cromwell_url', metavar='CROMWELL', type=url_regex_validator, nargs=1, help='Cromwell host')
    parser.add_argument('gcs_path', metavar='GCSPATH', type=gcs_path_regex_validator, nargs=1, help='GCS path to upload to')
    parser.add_argument('workflows', metavar='WORKFLOW', type=workflow_regex_validator, nargs='+', help='Workflows to process')

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()

    cromwell_url = args.cromwell_url[0]
    gcs_bucket, gcs_path = args.gcs_path[0]
    workflows = args.workflows

    credentials, project_id = google.auth.default()
    storage_client = storage.Client(credentials = credentials)
    papi_clients = PapiClients(credentials)

    logger.info(f'cromwell: {cromwell_url}')
    logger.info(f'gcs_bucket: {gcs_bucket}; gcs_path: {gcs_path}')
    logger.info(f'workflows: {workflows}')

    for workflow in workflows:
        process_workflow(cromwell_url, gcs_bucket, gcs_path, storage_client, papi_clients, workflow)

    logger.info('Extractor operation completed successfully.')
