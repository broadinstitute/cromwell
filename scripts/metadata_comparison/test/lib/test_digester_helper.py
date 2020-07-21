from google.cloud import storage
import logging
from metadata_comparison.lib.comparison_paths import ComparisonPath
from typing import AnyStr


VERSION_PAPI_V1 = 'PAPIv1'
VERSION_PAPI_V2 = 'PAPIv2_alpha1'


def download_metadata_from_gcs(bucket: storage.Bucket, local_sample_path: ComparisonPath) -> None:
    (local_sample_path / "operations").mkdir_p()

    prefix = str(local_sample_path)
    blobs = bucket.list_blobs(prefix=prefix)
    for blob in blobs:
        if not blob.name.endswith('/digest.json'):
            logging.info(f'Downloading blob: {blob.name}')
            blob.download_to_filename(blob.name)


def subdir_for_papi_version(papi_version: AnyStr) -> AnyStr:
    if papi_version == VERSION_PAPI_V1:
        path_element = 'PAPIv1'
    elif papi_version == VERSION_PAPI_V2:
        path_element = 'PAPIv2_alpha1/v1_style_machine_types'
    else:
        raise ValueError(f'Unrecognized PAPI version {papi_version}')
    return f'exome_germline_single_sample_v1.3/{path_element}'


def download_metadata_from_gcs_if_needed(
        sample_name: AnyStr, local_parent: ComparisonPath, bucket: storage.Bucket) -> None:
    """
    Copy down workflow and PAPI operations metadata from GCS if needed to test Local.
    """
    local_sample_path = local_parent / sample_name
    if not local_sample_path.exists():
        logging.info(f"Local sample directory '{local_sample_path}' does not exist, downloading from GCS.")
        download_metadata_from_gcs(bucket, local_sample_path)


def gcs_parent(subdir: AnyStr, gcs_comparison_path_by_subdir: dict) -> ComparisonPath:
    """
    GcsComparisonPaths are somewhat expensive to create so cache them.
    """
    if subdir not in gcs_comparison_path_by_subdir:
        path = ComparisonPath.create(f'gs://papi-performance-analysis/{subdir}')
        gcs_comparison_path_by_subdir[subdir] = path
    return gcs_comparison_path_by_subdir[subdir]
