#!/usr/bin/env python3

from google.cloud import storage
import logging

def set_log_verbosity(verbose: bool) -> None:
    if verbose:
        logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=logging.WARNING)


def quieten_chatty_imports() -> None:
    logging.getLogger('googleapiclient.discovery_cache').setLevel(logging.ERROR)
    logging.getLogger('googleapiclient.discovery').setLevel(logging.WARNING)


def upload_blob(bucket_name: str, source_file_contents: str, destination_blob_name: str, gcs_storage_client: storage.Client, logger: logging.Logger) -> None:
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"
    bucket = gcs_storage_client.bucket(bucket_name)
    logger.info(f'Uploading file content to gs://{bucket_name}/{destination_blob_name}...')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)
