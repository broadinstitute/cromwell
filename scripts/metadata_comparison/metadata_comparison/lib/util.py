import logging as log


def set_log_verbosity(verbose):
    level = log.INFO if verbose else log.WARNING
    log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=level)


def upload_blob(bucket_name, source_file_contents, destination_blob_name, gcs_storage_client, logger):
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"

    bucket = gcs_storage_client.bucket(bucket_name)

    logger.info(f'Uploading file content to gs://{bucket_name}/{destination_blob_name}...')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)


def quieten_chatty_imports():
    log.getLogger('googleapiclient.discovery_cache').setLevel(log.ERROR)
    log.getLogger('googleapiclient.discovery').setLevel(log.WARNING)


def ensure_slashed(path):
    return path + '/' if not path.endswith('/') else path


def build_papi_operation_mapping(json_metadata, call_fn):
    """Finds all instances of attempts with PAPI operations IDs in a workflow and
       invokes the supplied call_fn on each, potentially populating the returned dictionary."""
    # {
    #   "calls": {
    #     "workflow_name.task_name": [
    #       {
    #         "backend": "Papi"
    #         "jobId": "projects/broad-dsde-cromwell-dev/operations/5788303950667477684",
    papi_operation_mapping = {}

    def examine_calls(calls, path_so_far):
        for callname in calls:
            attempts = calls[callname]
            for attempt in attempts:
                operation_id = attempt.get('jobId')
                sub_workflow_metadata = attempt.get('subWorkflowMetadata')
                path = build_call_path(callname, path_so_far, attempt)
                if operation_id:
                    call_fn(papi_operation_mapping, operation_id, path, attempt)
                if sub_workflow_metadata:
                    examine_calls(sub_workflow_metadata.get('calls', {}), path)

    def build_call_path(callname, path_so_far, attempt):
        call_path = path_so_far.copy()

        # Remove confusing duplication in subworkflow call names
        deduplicated_callname = callname
        if len(path_so_far) > 0:
            this_call_components = callname.split('.')
            if len(this_call_components) > 1 and path_so_far[-1].endswith('.' + this_call_components[0]):
                deduplicated_callname = '.'.join(this_call_components[1:])

        call_path.append(deduplicated_callname)
        shard_index = attempt.get('shardIndex', -1)
        if shard_index != -1:
            call_path.append(f"shard_{shard_index:04d}")

        return call_path

    examine_calls(calls=json_metadata.get('calls', {}), path_so_far=[])

    return papi_operation_mapping
