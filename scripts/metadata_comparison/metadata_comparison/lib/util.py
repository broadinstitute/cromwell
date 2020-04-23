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


def ensure_slashed(string):
    return string + '/' if not string.endswith('/') else string
