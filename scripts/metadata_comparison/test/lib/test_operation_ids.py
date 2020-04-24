#!/usr/bin/env python3

import unittest
import json
from metadata_comparison.lib.operation_ids import *
from metadata_comparison.extractor import find_operation_ids_in_metadata
from test.lib.helper_functions import RESOURCES

class OperationIdTestMethods(unittest.TestCase):

    v1alpha2_ids = [
        'operations/EMj9o52aLhj78ZLxzunkiHcg0e2BmaAdKg9wcm9kdWN0aW9uUXVldWU',
        'operations/EMLel52aLhjHt6LH_fmuqOUBINHtgZmgHSoPcHJvZHVjdGlvblF1ZXVl',
        'operations/EJGZrp2aLhix7YG3sMaj-usBINHtgZmgHSoPcHJvZHVjdGlvblF1ZXVl',
        'operations/EM79o52aLhisyJifgbTuzb4BINHtgZmgHSoPcHJvZHVjdGlvblF1ZXVl'
    ]

    v2alpha1_ids = [
        'projects/broad-dsde-cromwell-dev/operations/4960504346170163809',
        'projects/broad-dsde-cromwell-dev/operations/1242302480525595574',
        'projects/broad-dsde-cromwell-dev/operations/11113854067012168443',
        'projects/broad-dsde-cromwell-dev/operations/14350975406210565808'
    ]

    v2beta_ids = [
        'projects/1005074806481/locations/us-central1/operations/14642412977689025366',
        'projects/1005074806481/locations/us-central1/operations/3128777609532869613',
        'projects/1005074806481/locations/us-central1/operations/18107476451113522273',
        'projects/1005074806481/locations/us-central1/operations/13032426715870634389'
    ]

    def test_operation_id_number_valid(self):
        v1alpha2_short_ids = [
            'EMj9o52aLhj78ZLxzunkiHcg0e2BmaAdKg9wcm9kdWN0aW9uUXVldWU',
            'EMLel52aLhjHt6LH_fmuqOUBINHtgZmgHSoPcHJvZHVjdGlvblF1ZXVl',
            'EJGZrp2aLhix7YG3sMaj-usBINHtgZmgHSoPcHJvZHVjdGlvblF1ZXVl',
            'EM79o52aLhisyJifgbTuzb4BINHtgZmgHSoPcHJvZHVjdGlvblF1ZXVl'
        ]
        v1_cases = zip(self.v1alpha2_ids, v1alpha2_short_ids)
        v2alpha1_short_ids = [
            '4960504346170163809',
            '1242302480525595574',
            '11113854067012168443',
            '14350975406210565808'
        ]
        v2alpha1_cases = zip(self.v2alpha1_ids, v2alpha1_short_ids)
        v2beta_short_ids = [
            '14642412977689025366',
            '3128777609532869613',
            '18107476451113522273',
            '13032426715870634389'
        ]
        v2beta_cases = zip(self.v2beta_ids, v2beta_short_ids)
        all_cases = list(v1_cases) + list(v2alpha1_cases) + list(v2beta_cases)

        for case in all_cases:
            with self.subTest(case=case):
                input = case[0]
                expectation = case[1]
                self.assertEqual(get_operation_id_number(input), expectation)

    def test_find_v1alpha2_operation_ids_in_metadata(self):
        metadata = json.loads((RESOURCES / 'forkjoin_metadata_papi_v1alpha2.json').read_text())
        self.assertEqual(find_operation_ids_in_metadata(metadata), self.v1alpha2_ids)


    def test_find_v2alpha1_operation_ids_in_metadata(self):
        metadata = json.loads((RESOURCES / 'forkjoin_metadata_papi_v2alpha1.json').read_text())
        self.assertEqual(find_operation_ids_in_metadata(metadata), self.v2alpha1_ids)


    def test_find_v2beta_operation_ids_in_metadata(self):
        metadata = json.loads((RESOURCES / 'forkjoin_metadata_papi_v2beta.json').read_text())
        self.assertEqual(find_operation_ids_in_metadata(metadata), self.v2beta_ids)


    def test_find_operation_ids_in_metadata_subworkflows(self):
        metadata = json.loads((RESOURCES / 'subworkflow_hello_world_metadata.json').read_text())
        self.assertEqual(find_operation_ids_in_metadata(metadata),
            ['projects/broad-dsde-cromwell-dev/operations/2244029211726316446'])


    def test_operation_id_to_api_version(self):
        for case in self.v1alpha2_ids:
            with self.subTest(case=case):
                self.assertEqual(operation_id_to_api_version(case), 'v1alpha2')
        for case in self.v2alpha1_ids:
            with self.subTest(case=case):
                self.assertEqual(operation_id_to_api_version(case), 'v2alpha1')
        for case in self.v2beta_ids:
            with self.subTest(case=case):
                self.assertEqual(operation_id_to_api_version(case), 'v2beta')


    def test_invalid_operation_id_has_no_api_version(self):
        case = "badstart/projects/broad-dsde-cromwell-dev/operations/4960504346170163809"
        with self.assertRaises(Exception) as context:
            operation_id_to_api_version(case)
        self.assertEqual(str(context.exception), 'Cannot deduce PAPI api version from unexpected operation ID format \'badstart/projects/broad-dsde-cromwell-dev/operations/4960504346170163809\'')



if __name__ == '__main__':
    unittest.main()
