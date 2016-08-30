# Introduction
Centaur is an integration testing suite for the [Cromwell](http://github.com/broadinstitute/cromwell) execution server.  It's purpose is to exercise the functionality of a specific deployment of Cromwell, to ensure that it is functioning properly 'in the wild'.  

# Prerequisites

Centaur expects to find a Cromwell server properly configured and running in server mode, listening on port 8000.  This can be configured by modifying the `cromwellUrl` parameter in `application.conf`.

# Running

There are two ways to invoke the intergration tests:
* `sbt test` - compiles and run via sbt directly, simple but also has the problem of running 2*cores tests in parallel which can overwhelm your Cromwell server if running in a development environment
* `run_tests_parallel.sh [THREADS]` - runs the same tests with an enforced parallelism limit.  Defaults to `3` if not specified

## Tags

All tests are tagged with their name and their TESTFORMAT, and also any custom tags specified in the `.test` file.
To run only those tests which have been tagged with a specified tag `tagFoo`:
```
sbt "test-only * -- -n tagFoo"
```

Or to instead exclude all tests which have been tagged with a specified tag `tagFoo`:
```
sbt "test-only * -- -l tagFoo"
```

# Adding custom tests

You can add your own tests to the test suite by adding `-Dcentaur.optionalTestPath=DIR` on your sbt invocation, 
e.g. `sbt -Dcentaur.optionalTestPath=/some/path/to/tests test`. The value of `DIR` is expected to be a directory
which contains one or more test case files.
 
The same result can be achieved more permanently by adding the custom directory into the application.conf file directly: 
```
centaur {
  optionalTestPath = "/some/path/to/tests"
}
```

# Defining test cases

Each test case file is a HOCON file with the following structure:
```
name: NAME  // Required: Name of the test
testFormat: TESTFORMAT // Required: One of WorkflowSuccessTest, WorkflowFailureTest
backends: [BACKENDNAME1, BACKENDNAME2, ...] // Optional list of backends. If supplied, this test will be ignored if these backends are not supported by the Cromwell server
basePath: /an/optional/field  // Optional, location for the files {} entries to be found relative to
tags: [ "any", "custom", "tags" ]  // Optional, a set of custom tags to apply to this test
ignore: false  // Optional, whether centaur will ignore this test when running

files {
  wdl: path/to/wdl  // Required: path to the WDL file to submit
  inputs: optional/path/to/inputs  // Optional, a path to an inputs JSON to include in the submission
  options: optional/path/to/options  // Optional, a path to an options JSON to include in the submission
}

// Optional, some metadata to verify on workflow completion:
metadata {
  fully.qualified.key.name1: VALUE1
  fully.qualified.key.name2: VALUE2
}
```

The tags are optional. If supplied they will allow people to turn on or off this test case by including or excluding tags when running (see above).

The `basePath` field is optional, but if supplied all paths will be resolved from that directory. If it is not supplied, all paths will be resolved from the directory the test case file is in.

The `testFormat` field can be one of the following, case insensitive:
* `workflowsuccess`: The workflow being supplied is expected to successfully complete
* `workflowfailure`: The workflow being supplied is expected to fail

The `metadata` is optional. If supplied, Centaur will retrieve the metadata from the successfully completed workflow and compare the values retrieved to those supplied. At the moment the only fields supported are strings, numbers and booleans.
For any metadata values which require workflow ID (i.e, file paths), use <\<UUID>> as a placeholder instead. For example:
* "calls.hello.hello.stdout": "gs://google-project/jes/root/wdl/<\<UUID>>/call-task/task-stdout.log"


