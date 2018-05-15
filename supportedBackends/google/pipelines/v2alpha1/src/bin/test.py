from google.cloud import storage
import sys


def create_empty_folder(bucket, path):
    gcs_client = storage.Client()
    bucket = gcs_client.get_bucket(bucket)
    blob = bucket.blob(path)
    blob.upload_from_string("", content_type="application/x-www-form-urlencoded;charset=UTF-8")


for line in sys.stdin.readlines():
    create_empty_folder("tjeandet-cromwell-execs", "test/" + line)