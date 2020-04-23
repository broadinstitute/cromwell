#!/usr/bin/env python3

import json
import unittest
from metadata_comparison.lib.argument_regex import *
from metadata_comparison.digester import find_succeeded_shards
from pathlib import Path

class DigesterTestMethods(unittest.TestCase):

    def read_resource(self, filename):
        path = f'test/resources/{filename}'
        with open(path, 'r') as file:
            data = file.read()
        return data


    def test_find_operation_ids_in_metadata(self):
        all_the_wfs = []
        samples = [
            'dev_D5327.NA12878',
            'dev_C862.NA19238',
            'dev_D5327.NA12891',
            'dev_D5327.NA12892',
            'dev_C1963.CHMI_CHMI3_Nex1'
        ]

        for path in Path('test/resources/exome_germline_single_sample_v1.3/PAPIv2_alpha1/v1_style_machine_types').rglob('workflow.json'):
            print(f'Hey I found path {path}')
            with open(path, 'r') as file:
                data = file.read()
                metadata = json.loads(data)
                find_succeeded_shards(metadata)
        #     {
        #         'projects/broad-dsde-cromwell-dev/operations/4960504346170163809': 'genomics',
        #         'projects/broad-dsde-cromwell-dev/operations/1242302480525595574': 'genomics',
        #         'projects/broad-dsde-cromwell-dev/operations/11113854067012168443': 'genomics',
        #         'projects/broad-dsde-cromwell-dev/operations/14350975406210565808': 'genomics'
        #     })


if __name__ == '__main__':
    unittest.main()
