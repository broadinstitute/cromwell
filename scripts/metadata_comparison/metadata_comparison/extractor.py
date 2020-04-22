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
from googleapiclient.discovery import build as google_client_build
import logging as log
from metadata_comparison.lib.argument_regex import *
from metadata_comparison.lib.operation_ids import *

logger = log.getLogger('metadata_comparison.extractor')

def set_log_verbosity(verbose):
    if verbose:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.INFO)
    else:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.WARNING)


def quieten_chatty_imports():
    log.getLogger('googleapiclient.discovery_cache').setLevel(log.ERROR)
    log.getLogger('googleapiclient.discovery').setLevel(log.WARNING)


def uploadLocalCheckout():
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


def fetch_raw_workflow_metadata(cromwell_url, workflow):
    url = f'{cromwell_url}/api/workflows/v1/{workflow}/metadata?expandSubWorkflows=true'
    logger.info(f'Fetching Cromwell metadata from {url}...')
    result = requests.get(url)
    return result.content, result.json()


def find_operation_ids_in_metadata(json_metadata):
    """Finds all instances of PAPI operations IDs in a workflow, and the API to call to retrieve metadata"""
    # {
    #   "calls": {
    #     "workflow_name.task_name": [
    #       {
    #         "backend": "Papi"
    #         "jobId": "projects/broad-dsde-cromwell-dev/operations/5788303950667477684",
    papi_operation_to_api_mapping = {}

    def find_operation_ids_in_calls(calls):
        for callname in calls:
            attempts = calls[callname]
            for attempt in attempts:
                operation_id = attempt.get('jobId')
                subWorkflowMetadata = attempt.get('subWorkflowMetadata')
                if operation_id:
                    api = 'lifesciences' if 'beta' in attempt.get('backend', '').lower() else 'genomics'
                    papi_operation_to_api_mapping[operation_id] = api
                if subWorkflowMetadata:
                    find_operation_ids_in_calls(subWorkflowMetadata.get('calls', {}))

    find_operation_ids_in_calls(json_metadata.get('calls', {}))

    return papi_operation_to_api_mapping


def read_papi_v2alpha1_operation_metadata(operation_id, api, genomics_v2alpha1_client):
    """Reads the operations metadata for a pipelines API v2alpha1 job ID. Returns a python dict"""
    logger.info(f'Reading PAPI v2alpha1 operation metadata for {operation_id}...')
    result = genomics_v2alpha1_client.projects().operations().get(name=operation_id).execute()
    return result


def upload_blob(bucket_name, source_file_contents, destination_blob_name, gcs_storage_client):
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"

    bucket = gcs_storage_client.bucket(bucket_name)

    logger.info(f'Uploading file content to gs://{bucket_name}/{destination_blob_name}...')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)


def upload_workflow_metadata_json(bucket_name, raw_workflow_metadata, workflow_gcs_base_path, gcs_storage_client):
    workflow_gcs_metadata_upload_path = f'{workflow_gcs_base_path}/metadata.json'
    upload_blob(bucket_name, raw_workflow_metadata, workflow_gcs_metadata_upload_path, gcs_storage_client)


def upload_operations_metadata_json(bucket_name, operation_id, operations_metadata, workflow_gcs_base_path, gcs_storage_client):
    """Uploads metadata to cloud storage, as json"""
    operation_upload_path = f'{workflow_gcs_base_path}/operations/{get_operation_id_number(operation_id)}.json'
    formatted_metadata = json.dumps(operations_metadata, indent=2)
    upload_blob(bucket_name, formatted_metadata, operation_upload_path, gcs_storage_client)


def process_workflow(cromwell_url, gcs_bucket, gcs_path, gcs_storage_client, genomics_v2alpha1_client, workflow):
    raw_metadata, json_metadata = fetch_raw_workflow_metadata(cromwell_url, workflow)
    workflow_gcs_base_path = f'{gcs_path}/{workflow}/extractor'

    operation_ids = find_operation_ids_in_metadata(json_metadata)
    for (id, api) in operation_ids.items():
        operation_metadata = read_papi_v2alpha1_operation_metadata(id, api, genomics_v2alpha1_client)
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
    genomics_v2alpha1_client = google_client_build('genomics', 'v2alpha1', credentials = credentials)

    logger.info(f'cromwell: {cromwell_url}')
    logger.info(f'gcs_bucket: {gcs_bucket}; gcs_path: {gcs_path}')
    logger.info(f'workflows: {workflows}')

    for workflow in workflows:
        process_workflow(cromwell_url, gcs_bucket, gcs_path, storage_client, genomics_v2alpha1_client, workflow)

    logger.info('Extractor operation completed successfully.')
