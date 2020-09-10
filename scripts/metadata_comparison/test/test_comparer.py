#!/usr/bin/env python3

from copy import deepcopy
import json
from metadata_comparison.lib.comparison_paths import ComparisonPath
from metadata_comparison.lib.operation_ids import JsonObject
from metadata_comparison.comparer import compare_jsons, csv_string_from_data
from test.lib.storage import RESOURCES
from typing import AnyStr, List, Tuple
import unittest


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
    def __compare_for_exome_germline_single_sample_list(json_1: JsonObject, json_2: JsonObject) -> List[List[AnyStr]]:
        return compare_jsons(json_1, json_2, "PAPIv1", "PAPIv2", [
            "ExomeGermlineSingleSample.AggregatedBamQC.",
            "ExomeGermlineSingleSample.BamToCram.", "ExomeGermlineSingleSample.BamToGvcf.VariantCalling.",
            "ExomeGermlineSingleSample.UnmappedBamToAlignedBam.", "ExomeGermlineSingleSample."])

    @staticmethod
    def __compare_for_exome_germline_single_sample(json_1: JsonObject, json_2: JsonObject) -> AnyStr:
        data = ComparerTestMethods.__compare_for_exome_germline_single_sample_list(json_1, json_2)
        actual = csv_string_from_data(data)
        return actual

    def test_compare_valid_jsons_no_disk_info(self) -> None:
        cases = [
            ("papiv1_version2_good.json", "papiv2_version2_good.json"),
        ]

        for case in cases:
            with self.subTest(case=case):
                json_1, json_2 = [json_from_path(p) for p in [case[0], case[1]]]
                actual = self.__compare_for_exome_germline_single_sample(json_1, json_2)
                expected = self.__read_resource("version3_comparer_on_version2_digests_good.csv")
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
                actual = self.__compare_for_exome_germline_single_sample_list(json_1, json_2)
                # Get the indexes for all the '% increase' columns.
                percent_column_indexes = [idx for idx, val in enumerate(actual[0]) if val == '% increase']

                # Skip the four header rows + blank line.
                for i in range(5, len(actual)):
                    for j in percent_column_indexes:
                        self.assertEqual(actual[i][j], "0.00%")

    first_digest = json.loads("""
        {
            "version": "0.0.2",
            "id": "e7933689-7891-44b1-a06c-0650fcd8f72c",
            "calls": {
                "foo": {
                    "delocalizationTimeSeconds": 6.280967,
                    "dockerImagePullTimeSeconds": 74.051255,
                    "localizationTimeSeconds": 359.455651,
                    "machineType": "n1-standard-1",
                    "otherTimeSeconds": 0.00,
                    "papiTotalTimeSeconds": 485.0,
                    "startupTimeSeconds": 33.132954,
                    "userCommandTimeSeconds": 11.998067
                }
            }
        }
        """)

    second_digest = json.loads("""
        {
            "version": "0.0.2",
            "workflowId": "c4002d69-06fb-4d9f-89f3-5b4acf333e7d",
            "calls": {
                "foo": {
                    "delocalizationTimeSeconds": 10.250919,
                    "dockerImagePullTimeSeconds": 89.463444,
                    "localizationTimeSeconds": 385.235749,
                    "machineType": "n1-standard-1",
                    "otherTimeSeconds": 18.655614000000128,
                    "papiTotalTimeSeconds": 548.142499,
                    "startupTimeSeconds": 34.831524,
                    "userCommandTimeSeconds": 9.705249
                }
            }
        }
        """)

    def setUp(self) -> None:
        self.updated_first_digest = deepcopy(self.first_digest)
        self.updated_second_digest = deepcopy(self.second_digest)
        super().setUp()

    def __run_negative_test(self, message_1: AnyStr, message_2: AnyStr):
        cases = [
            (self.updated_first_digest, self.second_digest, message_1),
            (self.first_digest, self.updated_second_digest, message_2),
        ]

        for case in cases:
            with self.subTest(case=case):
                first, second, message = case
                with self.assertRaises(ValueError) as cm:
                    self.__compare_for_exome_germline_single_sample_list(first, second)
                self.assertEqual(message, str(cm.exception))

    def test_different_versions(self):
        self.updated_first_digest['version'] = '0.0.1'
        self.updated_second_digest['version'] = '0.0.1'

        self.__run_negative_test('Inconsistent digest versions: First JSON digest is 0.0.1 but second is 0.0.2',
                                 'Inconsistent digest versions: First JSON digest is 0.0.2 but second is 0.0.1')

    def test_different_calls(self):
        self.updated_first_digest['calls']['bar'] = self.first_digest['calls']['foo']
        del self.updated_first_digest['calls']['foo']

        self.updated_second_digest['calls']['bar'] = self.second_digest['calls']['foo']
        del self.updated_second_digest['calls']['foo']

        message_1 = 'The specified digest files do not have the same call keys. ' + \
                    'These digests cannot be compared and probably are not from the same workflow and sample. ' + \
                    'In PAPIv1 but not PAPIv2: bar. In PAPIv1 but not PAPIv2: foo.'

        message_2 = 'The specified digest files do not have the same call keys. ' + \
                    'These digests cannot be compared and probably are not from the same workflow and sample. ' + \
                    'In PAPIv1 but not PAPIv2: foo. In PAPIv1 but not PAPIv2: bar.'

        self.__run_negative_test(message_1, message_2)

    def test_different_machine_types(self):
        self.updated_first_digest['calls']['foo']['machineType'] = 'n1-standard-2'
        self.updated_second_digest['calls']['foo']['machineType'] = 'n1-standard-2'

        message_1 = 'The specified digest files unexpectedly contain corresponding jobs that ran with different ' + \
                    'machine types: foo: n1-standard-2 vs n1-standard-1. ' + \
                    'Specify the --force argument to force comparison anyway.'

        message_2 = 'The specified digest files unexpectedly contain corresponding jobs that ran with different ' + \
                    'machine types: foo: n1-standard-1 vs n1-standard-2. ' + \
                    'Specify the --force argument to force comparison anyway.'

        self.__run_negative_test(message_1, message_2)

    def test_missing_required_key(self):
        del(self.updated_first_digest['calls']['foo']['machineType'])
        del(self.updated_second_digest['calls']['foo']['papiTotalTimeSeconds'])

        self.__run_negative_test(
             "In first digest JSON: call 'foo' missing required key 'machineType'",
             "In second digest JSON: call 'foo' missing required key 'papiTotalTimeSeconds'")


def json_from_path(string: AnyStr) -> ComparisonPath:
    path = ComparisonPath.create(str(RESOURCES) + '/comparer/' + string)
    return json.loads(path.read_text())


if __name__ == '__main__':
    unittest.main()
