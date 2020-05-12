#!/usr/bin/env python3

import os
import unittest
from typing import AnyStr, Callable
from metadata_comparison.lib.comparison_paths import ComparisonPath
from pathlib import Path
import logging
from metadata_comparison.digester import digest
from metadata_comparison.lib.logging import quieten_chatty_imports, set_log_verbosity


def read_resource(filename: AnyStr) -> AnyStr:
    path = Path('test/resources') / filename
    with open(path, 'r') as file:
        data = file.read()
    return data


class DigesterTestMethods(unittest.TestCase):
    set_log_verbosity(verbose=True)
    quieten_chatty_imports()

    def test_digestion(self) -> None:
        """
        This uses "real" metadata from the PAPI v2 performance spike to drive digester testing. The metadata is stored
        in GCS and copied down to the local machine if not already present from an earlier run. The digester can run
        against either local or GCS paths using `ComparisonPath`s. Local is nicer to iterate on than GCS since it
        runs so much more quickly; GCS testing can be turned off by setting the DIGESTER_TEST_LOCAL_ONLY environment
        variable.
        """
        subdir = 'exome_germline_single_sample_v1.3/PAPIv2_alpha1/v1_style_machine_types'
        local_parent = ComparisonPath.create(subdir)
        gcs_parent = ComparisonPath.create(f'gs://papi-performance-analysis/{subdir}')
        samples = {
            'dev_C1963.CHMI_CHMI3_Nex1': {
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
            # more samples if needed
            # 'dev_C862.NA19238',
            # 'dev_D5327.NA12878',
            # 'dev_D5327.NA12891',
            # 'dev_D5327.NA12892',
            # 'dev_RP-1535.NA17-308'
        }

        for sample_name in samples.keys():
            # Copy down workflow and PAPI operations metadata from GCS if needed to test Local
            local_sample_path = local_parent / sample_name
            if not local_sample_path.exists():
                logging.info(f"Local sample directory '{local_sample_path}' does not exist, downloading from GCS.")
                local_sample_path.mkdir()
                command = f"gsutil -m cp -r {gcs_parent}/{sample_name}/ {local_parent}"
                logging.info(f'Executing command: {command}')
                os.system(command)

            parents_to_test = [local_parent, gcs_parent]

            for parent in parents_to_test:
                description = parent.description()
                logging.info(f"Running digester test on {description} for sample '{sample_name}'")
                sample_path = parent / sample_name
                workflow_path = sample_path / 'workflow.json'
                operations_path = sample_path / 'operations'
                actual = digest(workflow_path, operations_path)

                expected = samples[sample_name]
                calls = actual.get('calls')

                actual_total = len(calls)
                self.assertEqual(actual_total, expected['total_jobs'])

                def more_than_x_attempts(attempts: int) -> Callable[[AnyStr], bool]:
                    def inner(call_name: AnyStr) -> bool:
                        return calls.get(call_name).get('attempt') > attempts
                    return inner

                for num_attempts in [1, 2, 3]:
                    actual = len(list(filter(more_than_x_attempts(num_attempts), calls)))
                    self.assertEqual(actual, expected[f'more_than_{num_attempts}_attempts'])

                def more_than_x_minutes_longer(minutes: int) -> Callable[[AnyStr], bool]:
                    def inner(call_name: AnyStr) -> bool:
                        return calls.get(call_name).get('cromwellAdditionalTotalTimeSeconds') > minutes * 60
                    return inner

                for minutes_longer in range(3, 9):
                    actual = len(list(filter(more_than_x_minutes_longer(minutes_longer), calls)))
                    self.assertEqual(actual, expected[f'cromwell_time_more_than_{minutes_longer}_minutes_longer_total'])


if __name__ == '__main__':
    unittest.main()
