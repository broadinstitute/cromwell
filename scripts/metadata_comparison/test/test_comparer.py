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

        print(f"expected dataframe:\n{expectedDf}")
        print(f"\nactual dataframe:\n{actualDf}")

        self.assertTrue(pandas.DataFrame.equals(expectedDf, actualDf))


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
