#!/usr/bin/env python3

import json
import os
import unittest
from typing import AnyStr
from metadata_comparison.lib.comparison_paths import ComparisonPath
from pathlib import Path
import logging
from metadata_comparison.lib.logging import quieten_chatty_imports, set_log_verbosity
from metadata_comparison.lib.operations_digesters import OperationDigester
import google.auth
from google.cloud import storage


class OperationsDigesterTestMethods(unittest.TestCase):
    set_log_verbosity(verbose=True)
    quieten_chatty_imports()

    def test_operations_digestion(self) -> None:
        """
        This uses "real" metadata from the PAPI v2 performance spike to drive operations digester testing.
        The metadata is stored in GCS and copied down to the local machine if not already present from an earlier run.
        Operations digesters can run against either local or GCS paths using `ComparisonPath`s. Since it is slow GCS
        testing is off by default, it can be turned on by setting the DIGESTER_TEST_GCS environment variable.
        """

        credentials, project_id = google.auth.default()
        storage_client = storage.Client(credentials=credentials)

        bucket_name = 'papi-performance-analysis'
        bucket = storage_client.get_bucket(bucket_name)

        # A cache of expensive-to-create GCS comparison paths.
        gcs_comparison_path_by_subdir = {}
        papi_versions = [VERSION_PAPI_V1, VERSION_PAPI_V2]

        for papi_version in papi_versions:
            subdir = subdir_for_papi_version(papi_version)
            local_parent = ComparisonPath.create(subdir)

            for sample_name in EXPECTATIONS.keys():
                download_metadata_from_gcs_if_needed(sample_name, local_parent, bucket)
                parents_to_test = [local_parent]
                # Skip slow GCS testing unless this environment variable is set.
                if os.environ.get('DIGESTER_TEST_GCS'):
                    parents_to_test.append(gcs_parent(subdir, gcs_comparison_path_by_subdir))

                for parent in parents_to_test:
                    description = parent.description()
                    logging.info(
                        f"Running operation digester on {description} sample '{sample_name}' backend {papi_version}")
                    sample_path = parent / sample_name

                    for operation in EXPECTATIONS.get(sample_name).get(papi_version).keys():
                        operations_path = sample_path / 'operations' / f'{operation}.json'
                        json_str = operations_path.read_text()
                        op_digester = OperationDigester.create(json.loads(json_str))
                        for key, value in EXPECTATIONS.get(sample_name).get(papi_version).get(operation).items():
                            method_to_call = getattr(op_digester, key)
                            self.assertEqual(method_to_call(), value, f'{key} was not {value}')


def read_resource(filename: AnyStr) -> AnyStr:
    path = Path('test/resources') / filename
    with open(path, 'r') as file:
        data = file.read()
    return data


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


VERSION_PAPI_V1 = 'PAPIv1'
VERSION_PAPI_V2 = 'PAPIv2_alpha1'

EXPECTATIONS = {
    'dev_C1963.CHMI_CHMI3_Nex1': {
        'PAPIv1': {
            'EODRwtmGLhiUmLjJhoTsspABIMu2sLePDSoPcHJvZHVjdGlvblF1ZXVl': {
            }
        },
        'PAPIv2_alpha1': {
            '9990846134018347343': {
                'total_time_seconds': 164.94303,
                'startup_time_seconds': 38.283515,
                'docker_image_pull_time_seconds': 63.978705,
                'localization_time_seconds': 23.193772,
                'user_command_time_seconds': 2.978336,
                'delocalization_time_seconds': 13.615169
            }
        }
    }
    # more samples if needed
    # 'dev_C862.NA19238',
    # 'dev_D5327.NA12878',
    # 'dev_D5327.NA12891',
    # 'dev_D5327.NA12892',
    # 'dev_RP-1535.NA17-308'
}


def gcs_parent(subdir: AnyStr, gcs_comparison_path_by_subdir: dict) -> ComparisonPath:
    """
    GcsComparisonPaths are somewhat expensive to create so cache them.
    """
    if subdir not in gcs_comparison_path_by_subdir:
        path = ComparisonPath.create(f'gs://papi-performance-analysis/{subdir}')
        gcs_comparison_path_by_subdir[subdir] = path
    return gcs_comparison_path_by_subdir[subdir]


if __name__ == '__main__':
    unittest.main()
