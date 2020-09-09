#!/usr/bin/env python3

import google.auth
from google.cloud import storage
import json
import logging
from metadata_comparison.lib.comparison_paths import ComparisonPath
from metadata_comparison.lib.logging import quieten_chatty_imports, set_log_verbosity
from metadata_comparison.lib.operations_digesters import Disk, DiskType, OperationDigester
import os
import unittest
from test.lib.test_digester_helper import download_metadata_from_gcs_if_needed,\
    gcs_parent, subdir_for_papi_version, VERSION_PAPI_V1, VERSION_PAPI_V2
from typing import Set


class OperationsDigesterTestMethods(unittest.TestCase):
    set_log_verbosity(verbose=True)
    quieten_chatty_imports()

    def test_operations_digestion(self) -> None:
        """
        This uses "real" metadata from the PAPI v2 performance spike to drive operations digester testing.
        The metadata is stored in GCS and copied down to the local machine if not already present from an earlier run.
        Operations digesters can run against either local or GCS paths using `ComparisonPath`s. Since GCS testing is
        slow it's turned off by default, it can be turned on by setting the DIGESTER_TEST_GCS environment variable.
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
                            if key == 'disks':
                                self.assertEquals(value, op_digester.disks())
                            else:
                                method_to_call = getattr(op_digester, key)
                                self.assertEqual(method_to_call(), value, f'{key} was not {value}')


EXPECTATIONS = {
    'dev_C1963.CHMI_CHMI3_Nex1': {
        'PAPIv1': {
            'EODRwtmGLhiUmLjJhoTsspABIMu2sLePDSoPcHJvZHVjdGlvblF1ZXVl': {
                'total_time_seconds': 466,
                'startup_time_seconds': 34.263103,
                'docker_image_pull_time_seconds': 88.822075,
                'localization_time_seconds': 47.902792,
                'user_command_time_seconds': 286.993021,
                'delocalization_time_seconds': 8.452244,
                'disks': {
                    'boot-disk': {
                        'name': 'boot-disk',
                        'sizeGb': 10,
                        'type': 'HDD'
                    },
                    'local-disk': {
                        'name': 'local-disk',
                        'sizeGb': 25,
                        'type': 'HDD'
                    }
                }
            }
        },
        'PAPIv2_alpha1': {
            '9990846134018347343': {
                'total_time_seconds': 164.94303,
                'startup_time_seconds': 38.283515,
                'docker_image_pull_time_seconds': 63.978705,
                'localization_time_seconds': 23.193772,
                'user_command_time_seconds': 2.978336,
                'delocalization_time_seconds': 13.615169,
                'disks': {
                    'boot-disk': {
                        'name': 'boot-disk',
                        'sizeGb': 11,
                        'type': 'HDD'
                    },
                    'local-disk': {
                        'name': 'local-disk',
                        'sizeGb': 10,
                        'type': 'HDD'
                    }
                }
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


if __name__ == '__main__':
    unittest.main()
