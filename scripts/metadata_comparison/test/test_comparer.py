#!/usr/bin/env python3

import json
import unittest
from typing import AnyStr, List, Tuple
from metadata_comparison.lib.comparison_paths import ComparisonPath
from metadata_comparison.lib.operation_ids import JsonObject
from metadata_comparison.comparer import compare_jsons, csv_string_from_data
from test.lib.storage import RESOURCES


class ComparerTestMethods(unittest.TestCase):

    @staticmethod
    def __read_resource(filename: str) -> AnyStr:
        return (ComparisonPath.create(str(RESOURCES)) / "comparer" / filename).read_text()

    @staticmethod
    def __read_test_json(filename: str) -> dict:
        return json.loads(ComparerTestMethods.__read_resource(filename))

    def __produce_workflow_id_and_json_tuples(self, workflow_id_and_filename_1: Tuple[str, str],
                                              workflow_id_and_filename_2: Tuple[str, str]) -> List[Tuple[str, dict]]:
        json1 = self.__read_test_json(workflow_id_and_filename_1[1])
        json2 = self.__read_test_json(workflow_id_and_filename_2[1])
        return [(workflow_id_and_filename_1[0], json1), (workflow_id_and_filename_2[0], json2)]

    @staticmethod
    def __compare_for_exome_germline_single_sample_list(json_1: JsonObject, json_2: JsonObject,
                                                        name_1: AnyStr, name_2: AnyStr) -> List[List[AnyStr]]:

        return compare_jsons(json_1, json_2, name_1, name_2, [
            "ExomeGermlineSingleSample.AggregatedBamQC.",
            "ExomeGermlineSingleSample.BamToCram.", "ExomeGermlineSingleSample.BamToGvcf.VariantCalling.",
            "ExomeGermlineSingleSample.UnmappedBamToAlignedBam.", "ExomeGermlineSingleSample."])

    @staticmethod
    def __compare_for_exome_germline_single_sample(json_1: JsonObject, json_2: JsonObject,
                                                   name_1: AnyStr, name_2: AnyStr) -> AnyStr:
        data = ComparerTestMethods.__compare_for_exome_germline_single_sample_list(json_1, json_2, name_1, name_2)
        actual = csv_string_from_data(data)
        return actual

    def test_compare_valid_jsons(self) -> None:
        cases = [
            ("papiv1_version2_good.json", "papiv2_version2_good.json"),
        ]

        for case in cases:
            with self.subTest(case=case):
                json_1, json_2 = [json_from_path(p) for p in [case[0], case[1]]]
                actual = self.__compare_for_exome_germline_single_sample(json_1, json_2, "PAPIv1", "PAPIv2")
                expected = self.__read_resource("good_comparison.csv")
                self.assertEqual(actual, expected)

    def test_compare_to_self(self) -> None:
        """ Sanity check comparing a digest to itself, which should produce 0.00% increases for everything. """
        cases = [
            ("papiv1_version2_good.json", "papiv1_version2_good.json"),
            ("papiv2_version2_good.json", "papiv2_version2_good.json"),
        ]

        for case in cases:
            with self.subTest(case=case):
                json_1, json_2 = [json_from_path(p) for p in [case[0], case[1]]]
                actual = self.__compare_for_exome_germline_single_sample_list(json_1, json_2, "PAPIv1", "PAPIv2")
                # Get the indexes for all the '% increase' columns.
                percent_column_indexes = [idx for idx, val in enumerate(actual[0]) if val == '% increase']
                print('foo')

                # Skip the four header rows + blank line.
                for i in range(5, len(actual)):
                    for j in percent_column_indexes:
                        self.assertEqual(actual[i][j], "0.00%")

#
#    def test_compare_invalid_jsons(self) -> None:
#        cases = [
#            ("performance_json_workflow_111.json", "performance_json_changed_key.json"),
#            ("performance_json_workflow_111.json", "performance_json_missing_key.json"),
#            ("performance_json_workflow_111.json", "performance_json_additional_key.json")
#        ]
#
#        for case in cases:
#            with self.subTest(case=case):
#                json_1, json_2 = [json_from_path(p) for p in [case[0], case[1]]]
#                with self.assertRaises(Exception) as context:
#                    compare_jsons(json_1, json_2)
#
#                self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


def json_from_path(string: AnyStr) -> ComparisonPath:
    path = ComparisonPath.create(str(RESOURCES) + '/comparer/' + string)
    return json.loads(path.read_text())


if __name__ == '__main__':
    unittest.main()
