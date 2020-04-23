#!/usr/bin/env python3

import unittest
import json
import pandas
from metadata_comparison.comparer import *

class ComparerTestMethods(unittest.TestCase):

    baseResourcePath = "test/resources/comparer"

    def readTestJson(self, filename):
        filePath = f"{self.baseResourcePath}/{filename}"
        with open(filePath,) as file:
            return filePath, json.load(file)


    def test_compare_valid_jsons(self):
        json1 = self.readTestJson("performance_json1.json")
        json2 = self.readTestJson("performance_json2.json")

        actualDf = compare_jsons(json1, json2)
        expectedDf = pandas.read_csv(f"{self.baseResourcePath}/valid_comparison_result.csv", index_col = 0)

        areEqual = pandas.DataFrame.equals(expectedDf, actualDf)
        if areEqual == False:
            # will print out dataframe having `true` in cells, which matching values and `false` otherwise
            print(expectedDf.eq(actualDf))

        self.assertTrue(areEqual)


    def test_compare_valid_differently_sorted_jsons(self):
        json1 = self.readTestJson("performance_json1.json")
        json2 = self.readTestJson("performance_json2_differently_sorted.json")

        # here we drop row names from dataframes `reset_index(drop = True)` in order to disregard file names during comparison
        actualDf = compare_jsons(json1, json2).reset_index(drop = True)
        expectedDf = pandas.read_csv(f"{self.baseResourcePath}/valid_comparison_result.csv", index_col = 0).reset_index(drop = True)

        areEqual = pandas.DataFrame.equals(expectedDf, actualDf)
        if areEqual == False:
            # will print out dataframe having `true` in cells, which matching values and `false` otherwise
            print(expectedDf.eq(actualDf))

        self.assertTrue(areEqual)


    def test_exception_on_renamed_key(self):
        json1 = self.readTestJson("performance_json1.json")
        json2 = self.readTestJson("performance_json_changed_key.json")
        with self.assertRaises(Exception) as context:
            compare_jsons(json1, json2)

        self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


    def test_exception_on_missing_key(self):
        json1 = self.readTestJson("performance_json1.json")
        json2 = self.readTestJson("performance_json_missing_key.json")
        with self.assertRaises(Exception) as context:
            compare_jsons(json1, json2)

        self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


    def test_exception_on_excessive_key(self):
        json1 = self.readTestJson("performance_json1.json")
        json2 = self.readTestJson("performance_json_excessive_key.json")
        with self.assertRaises(Exception) as context:
            compare_jsons(json1, json2)

        self.assertTrue("doesn't have matching subset of columns" in str(context.exception))


if __name__ == '__main__':
    unittest.main()
