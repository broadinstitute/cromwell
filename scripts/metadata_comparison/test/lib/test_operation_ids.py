#!/usr/bin/env python3

import unittest
import json
from metadata_comparison.lib.operation_ids import *
from metadata_comparison.extractor import find_operation_ids_in_metadata
from test.lib.helper_functions import read_resource

class OperationIdTestMethods(unittest.TestCase):

    def test_operation_id_number_valid(self):
        self.assertEqual(get_operation_id_number('projects/project_name/operations/01234567891011121314'), '01234567891011121314')


    def test_find_operation_ids_in_metadata(self):
        metadata = json.loads(read_resource('forkjoin_metadata.json'))
        self.assertEqual(find_operation_ids_in_metadata(metadata),
            [
                'projects/broad-dsde-cromwell-dev/operations/4960504346170163809',
                'projects/broad-dsde-cromwell-dev/operations/1242302480525595574',
                'projects/broad-dsde-cromwell-dev/operations/11113854067012168443',
                'projects/broad-dsde-cromwell-dev/operations/14350975406210565808'
            ])

    def test_find_operation_ids_in_metadata_subworkflows(self):
        metadata = json.loads(read_resource('subworkflow_hello_world_metadata.json'))
        self.assertEqual(find_operation_ids_in_metadata(metadata),
            ['projects/broad-dsde-cromwell-dev/operations/2244029211726316446'])


if __name__ == '__main__':
    unittest.main()
