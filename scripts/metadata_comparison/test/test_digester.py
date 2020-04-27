#!/usr/bin/env python3

import json
import unittest
from metadata_comparison.lib.argument_regex import *
from metadata_comparison.digester import *
from pathlib import Path


def read_resource(filename):
    path = f'test/resources/{filename}'
    with open(path, 'r') as file:
        data = file.read()
    return data


class DigesterTestMethods(unittest.TestCase):

    def test_digestion(self):
        for path in Path('test/resources/exome_germline_single_sample_v1.3/PAPIv2_alpha1/v1_style_machine_types').rglob('workflow.json'):
            parent = path.parent
            sample = parent.parts[-1]
            with open(path, 'r') as file:
                data = file.read()
                metadata = json.loads(data)
                actual = digest(metadata)
                with open(parent / 'digester' / Version / 'digest.json', 'r') as expected_file:
                    expected_data = expected_file.read()
                    expected = json.loads(expected_data)
                    self.assertEqual(actual, expected, f'oh noes {sample} has a mismatch')


if __name__ == '__main__':
    unittest.main()
