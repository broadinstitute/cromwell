#!/usr/bin/env python3

import json
import pandas
import unittest
from typing import AnyStr, List, Tuple
from pathlib import Path
from metadata_comparison.comparer import compare_jsons, ComparisonPath
from test.lib.storage import RESOURCES


class ComparerTestMethods(unittest.TestCase):

    valid_comparison_result_file = RESOURCES / Path("comparer/valid_comparison_result.csv")

    @staticmethod

    def __read_test_json(self, filename: str) -> dict:
        json_str = (RESOURCES / "comparer" / filename).read_text()
        return json.loads(json_str)

    def __produce_workflow_id_and_json_tuples(self, workflow_id_and_filename_1: Tuple[str, str], workflow_id_and_filename_2: Tuple[str, str]) -> List[Tuple[str, dict]]:
        json1 = self.__read_test_json(workflow_id_and_filename_1[1])
        json2 = self.__read_test_json(workflow_id_and_filename_2[1])
        return [(workflow_id_and_filename_1[0], json1), (workflow_id_and_filename_2[0], json2)]

    def test_compare_valid_jsons(self) -> None:
        cases = [
            ("performance_json_workflow_111.json", "performance_json_workflow_222.json"),
            ("performance_json_workflow_111.json", "performance_json_workflow_222_differently_sorted.json")
        ]

        for case in cases:
            with self.subTest(case=case):
                json_1, json_2 = [json_from_path(p) for p in [case[0], case[1]]]
                actual_df = compare_jsons(json_1, json_2)
                expected_df = pandas.read_csv(self.valid_comparison_result_file, index_col = 0)

                are_equal = pandas.DataFrame.equals(expected_df, actual_df)
                if not are_equal:
                    # will print out dataframe having `true` in cells, which matching values and `false` otherwise
                    print(expected_df.eq(actual_df))

                self.assertTrue(are_equal)

    def test_compare_invalid_jsons(self) -> None:
        cases = [
            ("performance_json_workflow_111.json", "performance_json_changed_key.json"),
            ("performance_json_workflow_111.json", "performance_json_missing_key.json"),
            ("performance_json_workflow_111.json", "performance_json_additional_key.json")
        ]

        for case in cases:
            with self.subTest(case=case):
                json_1, json_2 = [json_from_path(p) for p in [case[0], case[1]]]
                with self.assertRaises(Exception) as context:
                    compare_jsons(json_1, json_2)

                self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


def json_from_path(string: AnyStr) -> ComparisonPath:
    path = ComparisonPath.create(str(RESOURCES) + '/comparer/' + string)
    return json.loads(path.read_text())


if __name__ == '__main__':
    unittest.main()
