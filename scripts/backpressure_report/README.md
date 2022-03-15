# Backpressure Report

This `backpressure_report` Python project allows for measuring the amount of time Cromwell runner instances spend in a
high I/O state which triggers internal backpressure. While in this high I/O state the Cromwell runners will not hand out
job execution or restart check tokens, so job starts and restarts will be slowed until I/O returns to normal levels.

## Running

Installation:

Probably best done inside a [virtual environment](https://docs.python.org/3/library/venv.html)

```shell
pip install .
```

Usage:

```shell
python -m backpressure_report.main <files with JSON formatted Logs Explorer logs>
```

The program parses JSON-formatted Google Logs Explorer logs JSON looking for backpressure messages.
The Logs Explorer query should look like the following:

```
    resource.labels.container_name="cromwell1-runner-app"
    (jsonPayload.message=~"IoActor backpressure off" OR jsonPayload.message=~"Beginning IoActor backpressure")
```

Multiple input files may be required to capture logging output from an entire interval of interest since Google imposes
limits on the number of log entries that can be exported from a single query.

Output is a CSV file like:

```
Interval (1 hour),All pods,Pod 47z68,Pod 4hgd4,Pod 7svrs,Pod 9l2ld,Pod 9p9j4,Pod bj4vh,Pod d85vc,Pod gdp8x,Pod gth4r,Pod jkpbj,Pod jrgsx,Pod ltmvs,Pod mkdjt,Pod qt2bq,Pod th2p8,Pod thwz9,Pod xvcrk,Pod z7jfk
2022-01-01 05:00:00+00:00,62,20,0,0,0,42,0,0,0,0,0,0,0,0,0,0,0,0,0
2022-01-01 06:00:00+00:00,40,0,0,0,0,0,0,0,40,0,0,0,0,0,0,0,0,0,0
2022-01-01 07:00:00+00:00,20,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2022-01-01 08:00:00+00:00,40,20,0,0,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0
2022-01-01 09:00:00+00:00,110,0,0,0,0,70,0,0,40,0,0,0,0,0,0,0,0,0,0
...
```

The first column is the timestamp for the interval start, the second column is the sum of all backpressure durations from all runner
pods during that interval in seconds, and all subsequent columns are the backpressure durations for individual pods during the interval in seconds.

### Questions

- Q: Why not run the scripts directly, eg `python main.py`?
  - A: Running Python from this outer directory allows it to discover the `backpressure_report` 
  project, and thus allows imports across and between scripts.

## Unit tests

To run the Python unit tests from the top-level `backpressure_report` directory 
(ie the one containing this README.MD file), run:
```sh
python -m unittest discover -v
```

This will:
 - Find the `backpressure_report` project in that subdirectory.
   - And make it importable to other scripts.
 - Run the Python built-in unittest script, which will:
   - Discover the tests project in the `test` directory
   - Run them, verbosely.