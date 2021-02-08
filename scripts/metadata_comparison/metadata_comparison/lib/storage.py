#!/usr/bin/env python3

from google.cloud import storage
import logging
from typing import AnyStr


def upload_blob(bucket_name: str,
                source_file_contents: AnyStr,
                destination_blob_name: str,
                gcs_storage_client: storage.Client,
                logger: logging.Logger) -> None:
    """Uploads a file to the cloud"""
    # bucket_name = "your-bucket-name"
    # source_file_contents = "... some file contents..."
    # destination_blob_name = "storage/object/name"
    bucket = gcs_storage_client.bucket(bucket_name)
    logger.info(f'Uploading file content to gs://{bucket_name}/{destination_blob_name}...')
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_string(source_file_contents)
