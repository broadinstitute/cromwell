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
#   * pip3 install --upgrade gitpython
#
# Remember to login to create application default credentials before use:
#   % gcloud auth application-default login

import argparse
import json
import requests
from google.cloud import storage
import google.auth
from pathlib import Path
import git
import os
import zipfile
import logging
from metadata_comparison.lib.argument_regex import gcs_path_regex_validator, workflow_regex_validator
from metadata_comparison.lib.operation_ids import get_operation_id_number, visit_papi_operations, CallNameSequence, \
    JsonObject, OperationId
from metadata_comparison.lib.papi.papi_clients import PapiClients
from metadata_comparison.lib.storage import upload_blob
from typing import Any, AnyStr, List, Mapping, Sequence, Union
from metadata_comparison.lib.logging import quieten_chatty_imports, set_log_verbosity

logger = logging.getLogger('metadata_comparison.extractor')


def __create_snapshot_of_local_repo(repo: git.Repo, cromwell_snapshots_path: Union[Path, str]) -> Union[Path, str]:
    last_commit_hash = repo.head.commit.hexsha
    if not os.path.exists(cromwell_snapshots_path):
        os.makedirs(cromwell_snapshots_path)
    current_snapshot_path = cromwell_snapshots_path / last_commit_hash
    if not os.path.exists(current_snapshot_path):
        os.makedirs(current_snapshot_path)
        repo.clone(current_snapshot_path)
    return current_snapshot_path


def __create_zip_file(zip_file_path: Union[Path, str], current_snapshot_path: Union[Path, str]):
    with zipfile.ZipFile(zip_file_path, "a", allowZip64=False) as zip_file:
        for root, dirs, files in os.walk(current_snapshot_path):
            for file in files:
                zip_file.write(os.path.join(root, file))


def upload_local_checkout(cromwell_path: Path,
                          gcs_bucket: str,
                          gcs_path: str,
                          gcs_storage_client: storage.Client) -> None:
    cromwell_snapshots_path = cromwell_path.parent / "cromwell_snapshots"

    repo = git.Repo(cromwell_path)
    if repo.is_dirty():
        raise Exception("Unable to upload local checkout to GCS: repository is dirty - need to do check in first.")

    zip_file_name = f"cromwell_code.zip"
    zip_file_path = Path(cromwell_snapshots_path / zip_file_name)
    if not os.path.exists(zip_file_path):
        current_snapshot_path = __create_snapshot_of_local_repo(repo, cromwell_snapshots_path)
        __create_zip_file(zip_file_path, current_snapshot_path)

    upload_blob(gcs_bucket, zip_file_path.read_bytes(), f"{gcs_path}/{zip_file_name}", gcs_storage_client, logger)


def upload_local_config(config_path: Path, gcs_bucket: str, gcs_path: str, gcs_storage_client: storage.Client):
    configuration_file_name = "cromwell.conf"
    upload_blob(gcs_bucket, config_path.read_text(), f"{gcs_path}/{configuration_file_name}", gcs_storage_client, logger)


def fetch_raw_workflow_metadata(cromwell_url: str, workflow: str) -> (requests.Response, JsonObject):
    """Fetches workflow metadata for a workflow. Returns the raw response and the dict read from json"""
    url = f'{cromwell_url}/api/workflows/v1/{workflow}/metadata?expandSubWorkflows=true'
    logger.info(f'Fetching Cromwell metadata from {url}...')
    result = requests.get(url)
    return result.content, result.json()


def upload_workflow_metadata_json(bucket_name: str,
                                  raw_workflow_metadata: bytes,
                                  workflow_gcs_base_path: str,
                                  gcs_storage_client: storage.Client) -> None:
    workflow_gcs_metadata_upload_path = f'{workflow_gcs_base_path}/workflow.json'
    upload_blob(bucket_name, raw_workflow_metadata, workflow_gcs_metadata_upload_path, gcs_storage_client, logger)


def upload_operations_metadata_json(bucket_name: str,
                                    operation_id: str,
                                    operations_metadata: Mapping[str, Any],
                                    workflow_gcs_base_path: str,
                                    gcs_storage_client: storage.Client) -> None:
    """Uploads metadata to cloud storage, as json"""
    operation_upload_path = f'{workflow_gcs_base_path}/operations/{get_operation_id_number(operation_id)}.json'
    formatted_metadata = json.dumps(operations_metadata, indent=2)
    upload_blob(bucket_name, bytes(formatted_metadata, 'utf-8'), operation_upload_path, gcs_storage_client, logger)


def find_operation_ids_in_metadata(json_metadata: JsonObject) -> Sequence[AnyStr]:
    """Finds all instances of PAPI operations IDs in a workflow"""
    # Eg given:
    # {
    #   "calls": {
    #     "workflow_name.task_name": [
    #       {
    #         "jobId": "projects/broad-dsde-cromwell-dev/operations/01234567891011121314",
    # ...
    #
    # We want to extract "projects/broad-dsde-cromwell-dev/operations/01234567891011121314"
    def call_fn(acc: List[AnyStr],
                operation_id: OperationId,
                call_name_sequence: CallNameSequence,
                attempt: JsonObject) -> None:
        acc.append(operation_id)

    return visit_papi_operations(json_metadata, call_fn, initial_accumulator=[])


def process_workflow(cromwell_url: str,
                     gcs_bucket: str,
                     gcs_path: str,
                     gcs_storage_client: storage.Client,
                     papi_clients: PapiClients,
                     workflow: str) -> None:
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
    parser.add_argument('cromwell_url', metavar='CROMWELL', type=str, nargs=1,
                        help='Cromwell host')
    parser.add_argument('gcs_path', metavar='GCSPATH', type=gcs_path_regex_validator, nargs=1,
                        help='GCS path to upload to')
    parser.add_argument('workflows', metavar='WORKFLOW', type=workflow_regex_validator, nargs='+',
                        help='Workflows to process')
    parser.add_argument('cromwell_checkout_path', metavar='CROMWELLCHECKOUTPATH', type=Path,
                        help='Path to Cromwell git checkout used to run workflows')
    parser.add_argument('cromwell_config_path', metavar='CROMWELLCONFIGPATH', type=Path,
                        help='Path to Cromwell configuration file used to run workflows')

    args = parser.parse_args()
    set_log_verbosity(args.verbose)
    quieten_chatty_imports()

    cromwell_url = args.cromwell_url[0]
    gcs_bucket, gcs_path = args.gcs_path[0]
    workflows = args.workflows

    credentials, project_id = google.auth.default()
    storage_client = storage.Client(credentials=credentials)
    papi_clients = PapiClients(credentials)

    logger.info(f'cromwell: {cromwell_url}')
    logger.info(f'gcs_bucket: {gcs_bucket}; gcs_path: {gcs_path}')
    logger.info(f'workflows: {workflows}')

    for workflow in workflows:
        process_workflow(cromwell_url, gcs_bucket, gcs_path, storage_client, papi_clients, workflow)

    if args.cromwell_checkout_path:
        upload_local_checkout(args.cromwell_checkout_path, gcs_bucket, gcs_path, storage_client)
    if args.cromwell_config_path:
        upload_local_config(args.cromwell_config_path, gcs_bucket, gcs_path, storage_client)

    logger.info('Extractor operation completed successfully.')
