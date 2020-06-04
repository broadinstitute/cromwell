#!/usr/bin/env python3

import os
import unittest
from typing import AnyStr, Callable
from metadata_comparison.lib.comparison_paths import ComparisonPath
from pathlib import Path
import logging
from metadata_comparison.digester import digest
from metadata_comparison.lib.logging import quieten_chatty_imports, set_log_verbosity
from metadata_comparison.lib.operation_ids import JsonObject
import google.auth
from google.cloud import storage


class DigesterTestMethods(unittest.TestCase):
    set_log_verbosity(verbose=True)
    quieten_chatty_imports()

    def test_digestion(self) -> None:
        """
        This uses "real" metadata from the PAPI v2 performance spike to drive digester testing. The metadata is stored
        in GCS and copied down to the local machine if not already present from an earlier run. The digester can run
        against either local or GCS paths using `ComparisonPath`s. Local is nicer to iterate on than GCS since it
        runs so much more quickly. Since it is slow GCS testing is off by default, it can be turned on by setting the
        DIGESTER_TEST_GCS environment variable.
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
                        f"Running digester test on {description} for sample '{sample_name}' on backend {papi_version}")
                    sample_path = parent / sample_name
                    workflow_path = sample_path / 'workflow.json'
                    operations_path = sample_path / 'operations'
                    actual = digest(workflow_path, operations_path)

                    expected = EXPECTATIONS[sample_name][papi_version]
                    calls: JsonObject = actual.get('calls')

                    actual_total = len(calls)
                    self.assertEqual(actual_total, expected['total_jobs'])

                    for num_attempts in [1, 2, 3]:
                        actual_len = len(list(filter(more_than_x_attempts(calls, num_attempts), calls)))
                        self.assertEqual(actual_len, expected[f'more_than_{num_attempts}_attempts'])

                    for minutes_longer in range(3, 9):
                        actual_len = len(list(filter(more_than_x_minutes_longer(calls, minutes_longer), calls)))
                        expectation = expected[f'cromwell_time_more_than_{minutes_longer}_minutes_longer_total']
                        self.assertEqual(actual_len, expectation)

                    # Currently just a smoke test to assert not-completely-insane results for both v1 and v2 digesters.
                    self.assertTrue(all([calls[name].get('dockerImagePullSeconds') > 0 for name in calls]))


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


def download_metadata_from_gcs_if_needed(sample_name: AnyStr, local_parent: ComparisonPath, bucket: storage.Bucket) -> None:
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
            'total_jobs': 133,
            'more_than_1_attempts': 19,
            'more_than_2_attempts': 3,
            'more_than_3_attempts': 1,
            'cromwell_time_more_than_3_minutes_longer_total': 15,
            'cromwell_time_more_than_4_minutes_longer_total': 4,
            'cromwell_time_more_than_5_minutes_longer_total': 2,
            'cromwell_time_more_than_6_minutes_longer_total': 1,
            'cromwell_time_more_than_7_minutes_longer_total': 1,
            'cromwell_time_more_than_8_minutes_longer_total': 0,
        },
        'PAPIv2_alpha1': {
            'total_jobs': 133,
            'more_than_1_attempts': 12,
            'more_than_2_attempts': 1,
            'more_than_3_attempts': 0,
            'cromwell_time_more_than_3_minutes_longer_total': 21,
            'cromwell_time_more_than_4_minutes_longer_total': 7,
            'cromwell_time_more_than_5_minutes_longer_total': 4,
            'cromwell_time_more_than_6_minutes_longer_total': 2,
            'cromwell_time_more_than_7_minutes_longer_total': 1,
            'cromwell_time_more_than_8_minutes_longer_total': 0,
            # insert more intelligent assertions here
        }
    }
    # more samples if needed
    # 'dev_C862.NA19238',
    # 'dev_D5327.NA12878',
    # 'dev_D5327.NA12891',
    # 'dev_D5327.NA12892',
    # 'dev_RP-1535.NA17-308'
}


def more_than_x_attempts(calls: JsonObject, attempts: int) -> Callable[[AnyStr], bool]:
    """
    Return a function to filter the calls that had more than the specified number of attempts.
    """
    def inner(call_name: AnyStr) -> bool:
        return calls.get(call_name).get('attempt') > attempts

    return inner


def more_than_x_minutes_longer(calls: JsonObject, minutes: int) -> Callable[[AnyStr], bool]:
    """
    Return a function to filter the calls that ran for more than the specified number of minutes.
    """
    def inner(call_name: AnyStr) -> bool:
        return calls.get(call_name).get('cromwellAdditionalTotalTimeSeconds') > minutes * 60

    return inner


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
