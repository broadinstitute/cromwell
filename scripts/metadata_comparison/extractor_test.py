#!/usr/bin/env python3

import unittest
import json
import extractor

class ExtractorTestMethods(unittest.TestCase):

    def read_contents(self, file):
        with open(file, 'r') as file:
            data = file.read()
        return data



    def test_find_operation_ids_in_metadata(self):
        metadata = json.loads(self.read_contents('forkjoin_metadata.json'))
        self.assertEqual(extractor.find_operation_ids_in_metadata(metadata),
            {
                'projects/broad-dsde-cromwell-dev/operations/4960504346170163809': 'genomics',
                'projects/broad-dsde-cromwell-dev/operations/1242302480525595574': 'genomics',
                'projects/broad-dsde-cromwell-dev/operations/11113854067012168443': 'genomics',
                'projects/broad-dsde-cromwell-dev/operations/14350975406210565808': 'genomics'
            })

    def test_find_operation_ids_in_metadata_subworkflows(self):
        metadata = json.loads(self.read_contents('subworkflow_hello_world_metadata.json'))
        self.assertEqual(extractor.find_operation_ids_in_metadata(metadata),
            {'projects/broad-dsde-cromwell-dev/operations/2244029211726316446': 'genomics'})

if __name__ == '__main__':
    unittest.main()