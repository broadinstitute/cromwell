#!/usr/bin/env python3

import unittest
import json
from metadata_comparison.lib.argument_regex import *
from metadata_comparison.lib.operation_ids import *
from metadata_comparison.extractor import find_operation_ids_in_metadata

class ExtractorTestMethods(unittest.TestCase):

    def read_resource(self, filename):
        path = f'test/resources/{filename}'
        with open(path, 'r') as file:
            data = file.read()
        return data


    def test_cromwell_url_regex_valid(self):
         cases = [
             ('http://localhost', 'http://localhost'),
             ('http://localhost:8000', 'http://localhost:8000'),
             ('http://localhost:8000/', 'http://localhost:8000'),
             ('http://localhost:8000/api/workflows', 'http://localhost:8000'),
             ('http://localhost:8000/prefix/to/api/workflows', 'http://localhost:8000/prefix/to')
         ]

         for case in cases:
             with self.subTest(case=case):
                 input = case[0]
                 expectation = case[1]
                 self.assertEqual(url_regex_validator(input), expectation)


    def test_gcs_regex_valid(self):
        cases = [
            ('gs://bucket/path/to/directory', ('bucket', 'path/to/directory')),
            ('gs://bucket/path/to/directory', ('bucket', 'path/to/directory'))
        ]

        for case in cases:
            with self.subTest(case=case):
                input = case[0]
                expectation = case[1]
                self.assertEqual(gcs_path_regex_validator(input), expectation)


    def test_operation_id_number_valid(self):
        self.assertEqual(get_operation_id_number('projects/project_name/operations/01234567891011121314'), '01234567891011121314')


    def test_find_operation_ids_in_metadata(self):
        metadata = json.loads(self.read_resource('forkjoin_metadata.json'))
        self.assertEqual(find_operation_ids_in_metadata(metadata),
            {
                'projects/broad-dsde-cromwell-dev/operations/4960504346170163809': 'genomics',
                'projects/broad-dsde-cromwell-dev/operations/1242302480525595574': 'genomics',
                'projects/broad-dsde-cromwell-dev/operations/11113854067012168443': 'genomics',
                'projects/broad-dsde-cromwell-dev/operations/14350975406210565808': 'genomics'
            })

    def test_find_operation_ids_in_metadata_subworkflows(self):
        metadata = json.loads(self.read_resource('subworkflow_hello_world_metadata.json'))
        self.assertEqual(find_operation_ids_in_metadata(metadata),
            {'projects/broad-dsde-cromwell-dev/operations/2244029211726316446': 'genomics'})


if __name__ == '__main__':
    unittest.main()
