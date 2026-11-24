# DRS tests

## `drs_usa_jdr.inputs`

Used by `gcpbatch_drs_usa_jdr.test` and `gcpbatch_drs_usa_jdr_preresolve.test`.

DRSHub does not return service accounts for JDR paths. Therefore they shouldn't need to be localized using the
Cromwell custom DOS/DRS localizer.

However, the first file, `file1`, was generated before JDR started saving file names at the end of the `gsUri`.

Thus the JDR generated `gsUri` for `file1` is:

```
gs://broad-jade-dev-data-bucket/ca8edd48-e954-4c20-b911-b017fedffb67/585f3f19-985f-43b0-ab6a-79fa4c8310fc
```

With a `fileName` of `hello_jade.json`.

While the JDR generated `gsUri` for `file2` is:

```
gs://broad-jade-dev-data-bucket/e1941fb9-6537-4e1a-b70d-34352a3a7817/ad783b60-aeba-4055-8f7b-194880f37259/hello_jade_2.json
```

With a `fileName` of `hello_jade_2.json`.

The effect of this is that `file1` must always be downloaded by the Cromwell custom DOS/DRS localizer, which will be
sure to download the file to a path ending in the true file name `hello_jade.json`.

Meanwhile the `gsUri` for `file2` may be passed directly to specialized command line programs that can read the file
name and extension from the end of the GCS path.
