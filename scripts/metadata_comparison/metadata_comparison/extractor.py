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
import re
import json
import requests
from google.cloud import storage
import google.auth
from googleapiclient.discovery import build as google_client_build
import logging as log

logger = log.getLogger('metadata_comparison.extractor')

def set_log_verbosity(verbose):
    if verbose:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.INFO)
    else:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.WARNING)
        

def quieten_chatty_imports():
    log.getLogger('googleapiclient.discovery_cache').setLevel(log.ERROR)
    log.getLogger('googleapiclient.discovery').setLevel(log.WARNING)


def workflow_regex_validator(value):
    """Makes sure that a value is a valid Cromwell workflow ID then returns the workflow ID"""
    workflow_regex=re.compile('^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
    if not workflow_regex.match(value):
        msg = f'Invalid workflow ID {value}. Expected {workflow_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)
    else:
        return value


def url_regex_validator(value):
    """
    Validates then extract the root of the Cromwell URL from the various URL strings which might be provided.
    Deliberately flexible because it's tedious to remember which script requires which type of format.
    eg:
        'http://localhost' => 'http://localhost'
        'http://localhost:8000' => 'http://localhost:8000'
        'http://localhost:8000/' => 'http://localhost:8000'
        'http://localhost:8000/api/workflows/' => 'http://localhost:8000'
        'http://localhost:8000/custom/prefix/api/workflows/' => 'http://localhost:8000/custom/prefix'
    """
    url_regex = re.compile('(http(s?)://((?!/api).)*[^/])(/(api.*)?)?')
    m = url_regex.match(value)
    if m:
        return m.group(1)
    else:
        msg = f'Invalid Cromwell URL {value}. Expected {url_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)


def gcs_path_regex_validator(value):
    """
    Validates then extracts the bucket and object-path from a GS string. Returned as a pair.
    eg:
        'gs://bucket/path/to/directory/' -> ('bucket', 'path/to/directory')
    """
    gcs_regex = re.compile('^gs://([a-zA-Z0-9-]+)/(([a-zA-Z0-9-]+/)*[a-zA-Z0-9-]+)/?$')
    m = gcs_regex.match(value)
    if m:
        return m.group(1), m.group(2)
    else:
        msg = f'Invalid GCS path {value}. Expected {gcs_regex.pattern}'
        raise argparse.ArgumentTypeError(msg)


def get_operation_id_number(value):
    """
    Validates then extracts from PAPI operation IDs just the final number.
    eg:
        'projects/project_name/operations/01234567891011121314' -> '01234567891011121314'
    """

    operation_regex = re.compile('^.*/([^/]*)')
    m = operation_regex.search(value)
    if m:
        return m.group(1)
    else:
        msg = f'Unexpected operation ID {value}. Expected something like {operation_regex.pattern}'
        raise Exception(msg)


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

    calls = json_metadata.get('calls', {})
    for callname in calls:
        attempts = calls[callname]
        for attempt in attempts:
            operation_id = attempt.get('jobId')
            if operation_id:
                api = 'lifesciences' if 'beta' in attempt.get('backend', '').lower() else 'genomics'
                papi_operation_to_api_mapping[operation_id] = api
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
