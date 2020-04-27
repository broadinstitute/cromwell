import logging as log
from typing import Any, AnyStr, Callable, Dict, List, TypeVar


def set_log_verbosity(verbose) -> None:
    level = log.INFO if verbose else log.WARNING
    log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=level)


def upload_blob(bucket_name: AnyStr, source_file_contents: AnyStr, destination_blob_name: AnyStr,
                gcs_storage_client, logger):
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"

    bucket = gcs_storage_client.bucket(bucket_name)

    logger.info(f'Uploading file content to gs://{bucket_name}/{destination_blob_name}...')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)


def quieten_chatty_imports() -> None:
    log.getLogger('googleapiclient.discovery_cache').setLevel(log.ERROR)
    log.getLogger('googleapiclient.discovery').setLevel(log.WARNING)


K = TypeVar('K')  # Any type.
V = TypeVar('V')  # Any type.
# The operation mapping function chooses the type of both the keys and the values in the returned MappingDict.
MappingDict = Dict[K, V]

OperationMappingCallFunction = Callable[[MappingDict, AnyStr, List[AnyStr], Dict[AnyStr, Any]], None]


def build_papi_operation_mapping(json_metadata: AnyStr, call_fn: OperationMappingCallFunction) -> MappingDict:
    """Finds all instances of attempts with PAPI operations IDs in a workflow and
       invokes the supplied call_fn on each, potentially updating a returned dictionary."""
    # {
    #   "calls": {
    #     "workflow_name.task_name": [
    #       {
    #         "backend": "Papi"
    #         "jobId": "projects/broad-dsde-cromwell-dev/operations/5788303950667477684",
    papi_operation_mapping = {}

    def examine_calls(call_names: dict, path_so_far: List[AnyStr]) -> None:
        for call_name in call_names:
            attempts = call_names[call_name]
            for attempt in attempts:
                operation_id = attempt.get('jobId')
                sub_workflow_metadata = attempt.get('subWorkflowMetadata')
                path = build_call_path(call_name, path_so_far, attempt)
                if operation_id:
                    call_fn(papi_operation_mapping, operation_id, path, attempt)
                if sub_workflow_metadata:
                    examine_calls(sub_workflow_metadata.get('calls', {}), path)

    def build_call_path(call_name, path_so_far: List[AnyStr], attempt: dict) -> List[AnyStr]:
        call_path = path_so_far.copy()

        # Remove confusing duplication in subworkflow call names
        deduplicated_call_name = call_name
        if len(path_so_far) > 0:
            this_call_components = call_name.split('.')
            if len(this_call_components) > 1 and path_so_far[-1].endswith('.' + this_call_components[0]):
                deduplicated_call_name = '.'.join(this_call_components[1:])

        call_path.append(deduplicated_call_name)
        shard_index = attempt.get('shardIndex', -1)
        if shard_index != -1:
            call_path.append(f"shard_{shard_index:04d}")

        return call_path

    examine_calls(call_names=json_metadata.get('calls', {}), path_so_far=[])

    return papi_operation_mapping
