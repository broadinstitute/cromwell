#!/usr/bin/env python3

import unittest
from metadata_comparison.extractor import *

class ExtractorTestMethods(unittest.TestCase):

    def read_contents(self, file):
        with open(file, 'r') as file:
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


if __name__ == '__main__':
    unittest.main()