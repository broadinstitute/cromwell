#!/usr/bin/env python3

import json
import os
import unittest
from typing import AnyStr
from metadata_comparison.lib.comparison_paths import ComparisonPath
from pathlib import Path
import logging
from metadata_comparison.digester import Version, digest
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
        subdir = 'exome_germline_single_sample_v1.3/PAPIv2_alpha1/v1_style_machine_types'
        local_parent = ComparisonPath.create(subdir)
        gcs_parent = ComparisonPath.create(f'gs://papi-performance-analysis/{subdir}')
        samples = [
            'dev_C1963.CHMI_CHMI3_Nex1',
            # 'dev_C862.NA19238',
            # 'dev_D5327.NA12878',
            # 'dev_D5327.NA12891',
            # 'dev_D5327.NA12892',
            # 'dev_RP-1535.NA17-308'
        ]

        for sample in samples:
            # Copy down workflow and PAPI operations metadata from GCS if needed to test Local
            local_sample_path = local_parent / sample
            if not local_sample_path.exists():
                logging.info(f'Local sample directory does not exist, downloading from GCS.')
                local_sample_path.mkdir()
                command = f"gsutil -m cp -r {gcs_parent}/{sample}/ {local_parent}"
                logging.info(f'Executing command: {command}')
                os.system(command)

            for description, parent in [('local filesystem', local_parent), ('GCS', gcs_parent)]:
                logging.info(f"Running digester test on {description} for sample '{sample}'")
                sample_path = parent / sample
                workflow_path = sample_path / 'workflow.json'
                operations_path = sample_path / 'operations'
                actual = digest(workflow_path, operations_path)

                expected_file = sample_path / 'digests' / Version / 'digest.json'
                expected_data = expected_file.read_text()
                expected = json.loads(expected_data)

                self.assertEqual(actual, expected, f'oh noes {sample} has a mismatch on {description}')


if __name__ == '__main__':
    unittest.main()
