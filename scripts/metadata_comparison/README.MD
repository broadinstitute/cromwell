# Metadata Comparison Scripts

This `metadata_comparison` python project provides tools to compare workflows run
in different Cromwell environments to compare overall cost and performance. 

## Running a script

Choose a script to run. For this example we will use the `extractor`.

From this top-level directory `metadata_comparison` directory (ie the one 
containing this README.MD file), run:

```sh
# python3 -m metadata_comparison.extractor <ARGS>
```

### Questions

- Q: Why not run the scripts directly, eg `python3 extractor.py`?
  - A: Running python from this outer directory allows it to discover the `metadata_comparison` 
  project, and thus allows imports across and between scripts.

## Unit tests

To run the python unit tests from the top-level `metadata_comparison` directory 
(ie the one containing this README.MD file), run:
```sh
# python3 -m unittest discover -v
```

This will:
 - Find the `metadata_comparison` project in that subdirectory.
   - And make it importable to other scripts.
 - Run the python built-in unittest script, which will:
   - Discover the tests project in the `test` directory
   - Run them, verbosely.