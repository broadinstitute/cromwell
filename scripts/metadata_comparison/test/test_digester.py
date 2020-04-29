#!/usr/bin/env python3

import unittest
from metadata_comparison.digester import *
from pathlib import Path


def read_resource(filename: AnyStr) -> AnyStr:
    path = Path('test/resources') / filename
    with open(path, 'r') as file:
        data = file.read()
    return data


class DigesterTestMethods(unittest.TestCase):

    def test_digestion(self) -> None:
        directory = Path('test/resources/exome_germline_single_sample_v1.3/PAPIv2_alpha1/v1_style_machine_types')
        for workflow_file in directory.rglob('workflow.json'):
            parent = workflow_file.parent
            sample = parent.parts[-1]

            workflow_path = ComparisonPath.create(workflow_file)
            operations_path = ComparisonPath.create(parent) / "operations"
            actual = digest(workflow_path, operations_path)

            expected_file = parent / 'digests' / Version / 'digest.json'
            expected_data = expected_file.read_bytes()
            expected = json.loads(expected_data)

            self.assertEqual(actual, expected, f'oh noes {sample} has a mismatch')


if __name__ == '__main__':
    unittest.main()
