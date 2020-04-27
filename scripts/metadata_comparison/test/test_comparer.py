#!/usr/bin/env python3

import unittest
import json
import pandas
from typing import List, Tuple
from pathlib import Path
from metadata_comparison.comparer import compare_jsons
from test.lib.helper_functions import RESOURCES

class ComparerTestMethods(unittest.TestCase):

    valid_comparison_result_file = RESOURCES / Path("comparer/valid_comparison_result.csv")

    def __read_test_json(self, filename: str) -> dict:
        filePath = RESOURCES / "comparer" / filename
        with open(filePath,) as file:
            return json.load(file)


    def __produce_workflow_id_and_json_tuples(self, workflow_id_and_filename_1: Tuple[str, str], workflow_id_and_filename_2: Tuple[str, str]) -> List[Tuple[str, dict]]:
        json1 = self.__read_test_json(workflow_id_and_filename_1[1])
        json2 = self.__read_test_json(workflow_id_and_filename_2[1])
        return [(workflow_id_and_filename_1[0], json1), (workflow_id_and_filename_2[0], json2)]


    def test_compare_valid_jsons(self) -> None:
        actual_workflow_ids_and_jsons = self.__produce_workflow_id_and_json_tuples((111, "performance_json_workflow_111.json"), (222, "performance_json_workflow_222.json"))

        actual_df = compare_jsons(actual_workflow_ids_and_jsons)
        expected_df = pandas.read_csv(self.valid_comparison_result_file, index_col = 0)

        are_equal = pandas.DataFrame.equals(expected_df, actual_df)
        if are_equal == False:
            # will print out dataframe having `true` in cells, which matching values and `false` otherwise
            print(expected_df.eq(actual_df))

        self.assertTrue(are_equal)


    def test_compare_valid_differently_sorted_jsons(self) -> None:
        actual_workflow_ids_and_jsons = self.__produce_workflow_id_and_json_tuples((111, "performance_json_workflow_111.json"), (222, "performance_json_workflow_222_differently_sorted.json"))

        # here we drop row names from dataframes `reset_index(drop = True)` in order to disregard file names during comparison
        actual_df = compare_jsons(actual_workflow_ids_and_jsons).reset_index(drop = True)
        expected_df = pandas.read_csv(self.valid_comparison_result_file, index_col = 0).reset_index(drop = True)

        are_equal = pandas.DataFrame.equals(expected_df, actual_df)
        if are_equal == False:
            # will print out dataframe having `true` in cells, which matching values and `false` otherwise
            print(expected_df.eq(actual_df))

        self.assertTrue(are_equal)


    def test_exception_on_renamed_key(self) -> None:
        actual_workflow_ids_and_jsons = self.__produce_workflow_id_and_json_tuples((111, "performance_json_workflow_111.json"), (222, "performance_json_changed_key.json"))
        with self.assertRaises(Exception) as context:
            compare_jsons(actual_workflow_ids_and_jsons)

        self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


    def test_exception_on_missing_key(self) -> None:
        actual_workflow_ids_and_jsons = self.__produce_workflow_id_and_json_tuples((111, "performance_json_workflow_111.json"), (222, "performance_json_missing_key.json"))
        with self.assertRaises(Exception) as context:
            compare_jsons(actual_workflow_ids_and_jsons)

        self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


    def test_exception_on_excessive_key(self) -> None:
        actual_workflow_ids_and_jsons = self.__produce_workflow_id_and_json_tuples((111, "performance_json_workflow_111.json"), (222, "performance_json_excessive_key.json"))
        with self.assertRaises(Exception) as context:
            compare_jsons(actual_workflow_ids_and_jsons)

        self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


if __name__ == '__main__':
    unittest.main()
