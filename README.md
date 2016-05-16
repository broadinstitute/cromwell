# Introduction
Centaur is an integration testing suite for the [Cromwell](http://github.com/broadinstitute/cromwell) execution server.  It's purpose is to exercise the functionality of a specific deployment of Cromwell, to ensure that it is functioning properly 'in the wild'.  

# Prerequisites

Centaur expects to find a Cromwell server properly configured and running in server mode, listening on port 8000.  This can be configured by modifying the `cromwellUrl` parameter in `application.conf`.

# Running

There are two ways to invoke the intergration tests:
* `sbt test` - compiles and run via sbt directly, simple but also has the problem of running 2*cores tests in parallel which can overwhelm your Cromwell server if running in a development environment
* `run_tests_parallel.sh [THREADS]` - runs the same tests with an enforced parallelism limit.  Defaults to `3` if not specified

# Adding custom tests

You can add your own tests to the test suite by adding `-Dcentaur.optionalTestPath=DIR` on your sbt invocation, 
e.g. `sbt -Dcentaur.optionalTestPath=/some/path/to/tests test`. The value of `DIR` is expected to be a directory
which contains one or more test case files. Each test case file is a HOCON file with the following structure:
```
testName: TESTNAME
testFormat: TESTFORMAT
basePath: /an/optional/field

files {
  wdl: path/to/wdl
  inputs: optional/path/to/inputs
  options: optional/path/to/options
  metadata: optional/path/to/metadata
}
```

The `basePath` field is optional, but if supplied all paths will be resolved from that directory. If it is not supplied, all paths will be resolved from the directory the test case file is in.

The `testFormat` field can be one of the following, case insensitive:
* `workflowsuccess`: The workflow being supplied is expected to successfully complete
* `workflowfailure`: The workflow being supplied is expected to fail 
* `submissionfailure`: The workflow being supplied is expected to fail on workflow submission

# FAQs

##### I'm seeing a bunch of timeout errors in my tests
This is because you are trying to drive your poor Cromwell too hard!  Typically this happens when you are running everything on your local machine and are using a high level of parallelism (e.g. running via `sbt test`).  Reduce your parallelism by using `run_tests_parallel.sh` instead and/or beef up your configuration

![Image of a Centaur](https://upload.wikimedia.org/wikipedia/commons/3/34/Centaur_IV_tank_of_'H'_Troop,_2nd_Battery,_Royal_Marine_Armoured_Support_Group,_13_June_1944._B5457.jpg)

