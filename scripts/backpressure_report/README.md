# Backpressure Report

This `backpressure_report` Python project allows for measuring the amount of time Cromwell runner instances spend in a
high I/O state which triggers internal backpressure. While in this high I/O state the Cromwell runners will not hand out
job execution or restart check tokens, so job starts and restarts will be slowed until I/O returns to normal levels.

## Running

Usage:

```shell
python -m backpressure_report.main <files with JSON formatted Log Explorer logs>
```

The program parses JSON-formatted Google Logs Explorer logs JSON looking for backpressure messages.
The Logs Explorer query should look like the following:

```
    resource.labels.container_name="cromwell1-runner-app"
    (jsonPayload.message=~"IoActor backpressure off" OR jsonPayload.message=~"Beginning IoActor backpressure")
```

Multiple input files may be required to capture logging output from an entire interval of interest as Google imposes
limits on the number of log entries that can be exported from a single query.

### Questions

- Q: Why not run the scripts directly, eg `python main.py`?
  - A: Running Python from this outer directory allows it to discover the `backpressure_report` 
  project, and thus allows imports across and between scripts.

## Unit tests

To run the Python unit tests from the top-level `backpressure_report` directory 
(ie the one containing this README.MD file), run:
```sh
# python3 -m unittest discover -v
```

This will:
 - Find the `backpressure_report` project in that subdirectory.
   - And make it importable to other scripts.
 - Run the Python built-in unittest script, which will:
   - Discover the tests project in the `test` directory
   - Run them, verbosely.