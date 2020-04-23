#!/usr/bin/env python3
#
# comparer.py
#
# Purpose: Compare performance metadata JSON files produced by Digester and produce result in CSV format
#
# Usage: python3 comparer.py [-h] [-v] [--json_paths JSONPATH [JSONPATH ...]] [--output_path OUTPUTPATH]
#
# Python Prereqs (at least, the ones which I needed to manually install... YMMV):
#
#   * pip3 install --upgrade pandas
#

import argparse
import json
import pandas
from pathlib import Path
import logging as log

logger = log.getLogger('metadata_comparison.comparer')

def set_log_verbosity(verbose):
    if verbose:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.INFO)
    else:
        log.basicConfig(format='[%(asctime)s] [%(name)s] %(message)s', level=log.WARNING)


def read_json_files(*paths):
    """
    Produces list of tuples: [[path, json], [path, json], ...]
    """
    jsons = []
    for path in paths:
        jsons.append([path, json.load(open(path,))])

    return jsons

def compare_jsons(*pathsAndJsons):
    """
    Uses pandas library to convert JSONs into into dataframes, and concatenate those dataframes into a single one.
    Performs sanity check, producing exception, if at least one of the JSONs has doesn't have matching subset of keys.
    """
    result = pandas.DataFrame()
    lastCols = []
    for pathAndJson in pathsAndJsons:
        df = pandas.json_normalize(pathAndJson[1])
        cols = [c for c in df.columns if c[len(c)-8:] != '.attempt' or c == 'version']

        if lastCols and lastCols != cols:
            raise Exception(f"JSON data at {pathAndJson[0]} doesn't have matching subsets of columns. Expected: {lastCols} but got {cols}")

        lastCols = cols
        df.index = [pathAndJson[0]]
        result = pandas.concat([result, df[cols]])

    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Compare performance metadata JSONs and produce CSV result')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--json_paths', metavar='JSONPATH', type=Path, nargs='+', help='Paths to JSON files')
    parser.add_argument('--output_path', metavar='OUTPUTPATH', type=Path, nargs=1, help='Path for output CSV file')

    logger.info("Starting Comparer operation.")

    args = parser.parse_args()
    set_log_verbosity(args.verbose)

    pathsAndJsons = read_json_files(*args.json_paths)
    comparisonResultDf = compare_jsons(*pathsAndJsons)
    comparisonResultDf.to_csv(args.output_path[0])

    logger.info('Comparer operation completed successfully.')
