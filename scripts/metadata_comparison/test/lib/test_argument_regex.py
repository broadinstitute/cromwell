#!/usr/bin/env python3

import unittest
from metadata_comparison.lib.argument_regex import *
from test.lib.helper_functions import *

class ArgumentRegexTestMethods(unittest.TestCase):

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


if __name__ == '__main__':
    unittest.main()
