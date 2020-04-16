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
#   * pip3 install --upgrade google-api-python-client
#
# Remember to login to create application default credentials before use:
#   % gcloud auth application-default login

import argparse
import re
import json
import requests
from google.cloud import storage
import google.auth
from apiclient.discovery import build

def workflow_regex_validator(value):
    workflow_regex=re.compile('^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
    if not workflow_regex.match(value):
        msg = 'Invalid workflow ID {inv_id}. Expected {regex}'.format(inv_id=value, regex=workflow_regex.pattern)
        raise argparse.ArgumentTypeError(msg)
    else:
        return value


def url_regex_validator(value):
    url_regex = re.compile('(http(s?)://.*[^/])(/(api.*)?)?')
    m = url_regex.search(value)
    if m:
        return m.group(1)
    else:
        msg = 'Invalid Cromwell URL {inv_url}. Expected {regex}'.format(inv_url=value, regex=url_regex.pattern)
        raise argparse.ArgumentTypeError(msg)


def gcs_path_regex_validator(value):
    gcs_regex = re.compile('^gs://([a-zA-Z0-9-]+)/(([a-zA-Z0-9-]+/)*[a-zA-Z0-9-]+)/?$')
    m = gcs_regex.search(value)
    if m:
        return m.group(1), m.group(2)
    else:
        msg = 'Invalid GCS path {inv_gcs}. Expected {regex}'.format(inv_gcs=value, regex=gcs_regex.pattern)
        raise argparse.ArgumentTypeError(msg)


def get_operation_id_number(value):
    """Extracts from operation IDs like 'projects/project_name/operations/01234567891011121314' just the final
    number - eg '01234567891011121314'"""

    operation_regex = re.compile('^.*/([^/]*)')
    m = operation_regex.search(value)
    if m:
        return m.group(1)
    else:
        msg = 'Unexpected operation ID {inv_url}. Expected something like {regex}'.format(inv_url=value, regex=url_regex.pattern)
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
    url = '{url}/api/workflows/v1/{wfid}/metadata?expandSubWorkflows=true'.format(url=cromwell_url, wfid=workflow)
    print('Fetching Cromwell metadata from...', url, end='...')
    result = requests.get(url)
    print(' fetched!')
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
            jobId = attempt.get('jobId')
            if jobId:
                api = 'lifesciences' if 'beta' in attempt.get('backend', '').lower() else 'genomics'
                papi_operation_to_api_mapping[jobId] = api
    return papi_operation_to_api_mapping


def read_genomics_papiv1_operations_metadata(jobId, api, genomics_client):
    """Reads the operations metadata for a pipelines API v1 job."""
    return genomics_client.operations().get(name=jobId).execute()


def pretty_print_json(raw_json):
    return json.dumps(json.loads(raw_json), indent=4, separators=(',', ': '))


def read_operations_metadata(jobId, api, genomics_client):
    """Reads the operations metadata for a pipelines API v2alpha1 job ID. Returns a python dict"""
    # TODO: Reroute appropriately for different API versions
    print('Reading operation metadata for {jobId}...'.format(jobId=jobId), end='')
    result = genomics_client.projects().operations().get(name=jobId).execute()
    print(' done!')
    return result


def upload_blob(bucket_name, source_file_contents, destination_blob_name, gcs_storage_client):
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"

    bucket = gcs_storage_client.bucket(bucket_name)

    print('Uploading file content to gs://{bucket}/{dest}...'.format(bucket=bucket_name, dest=destination_blob_name), end='')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)
    print(" done!")


def upload_workflow_metadata_json(bucket_name, raw_workflow_metadata, workflow_gcs_base_path, gcs_storage_client):
    workflow_gcs_metadata_upload_path = '{workflow_gcs_base_path}/metadata.json'.format(workflow_gcs_base_path=workflow_gcs_base_path)
    upload_blob(bucket_name, raw_workflow_metadata, workflow_gcs_metadata_upload_path, gcs_storage_client)


def upload_operations_metadata_json(bucket_name, operation_id, operations_metadata, workflow_gcs_base_path, gcs_storage_client):
    """Uploads metadata to cloud storage, as json"""
    operation_upload_path = '{workflow_gcs_base_path}/operations/{miniid}.json'.format(workflow_gcs_base_path=workflow_gcs_base_path, miniid=get_operation_id_number(operation_id))
    formatted_metadata = json.dumps(operations_metadata, indent=2)
    upload_blob(bucket_name, formatted_metadata, operation_upload_path, gcs_storage_client)


def process_workflow(cromwell_url, gcs_bucket, gcs_path, gcs_storage_client, genomics_client, workflow):
    raw_metadata, json_metadata = fetch_raw_workflow_metadata(cromwell_url, workflow)
    workflow_gcs_base_path = '{root_gcs_path}/{wfid}/extractor'.format(root_gcs_path=gcs_path, wfid=workflow)

    operation_ids = find_operation_ids_in_metadata(json_metadata)
    operation_metadata = { id: read_operations_metadata(id, api, genomics_client) for (id, api) in operation_ids.items() }
    for id in operation_metadata:
        upload_operations_metadata_json(gcs_bucket, id, operation_metadata[id], workflow_gcs_base_path, gcs_storage_client)
    upload_workflow_metadata_json(gcs_bucket, raw_metadata, workflow_gcs_base_path, gcs_storage_client)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract metadata and operation details for workflows and upload to GCS')
    parser.add_argument('cromwell_url', metavar='CROMWELL', type=url_regex_validator, nargs=1, help='Cromwell host')
    parser.add_argument('gcs_path', metavar='GCSPATH', type=gcs_path_regex_validator, nargs=1, help='GCS path to upload to')
    parser.add_argument('workflows', metavar='WORKFLOW', type=workflow_regex_validator, nargs='+', help='Workflows to process')

    args = parser.parse_args()
    cromwell_url = args.cromwell_url[0]
    gcs_bucket, gcs_path = args.gcs_path[0]
    workflows = args.workflows

    credentials, project_id = google.auth.default()
    storage_client = storage.Client(credentials = credentials)
    genomics_v1_client = build('genomics', 'v1alpha2', credentials = credentials)
    genomics_client = build('genomics', 'v2alpha1', credentials = credentials)

    print('cromwell:', cromwell_url)
    print('gcs_bucket', gcs_bucket, 'gcs_path', gcs_path)
    print('workflows:', workflows)

    for workflow in workflows:
        process_workflow(cromwell_url, gcs_bucket, gcs_path, storage_client, genomics_client, workflow)
