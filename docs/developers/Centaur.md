Centaur is an integration testing suite for the [Cromwell](http://github.com/broadinstitute/cromwell) execution engine.  Its purpose is to exercise the functionality of a specific deployment of Cromwell, to ensure that it is functioning properly 'in the wild'.  

## Prerequisites

Centaur expects to find a Cromwell server properly configured and running in server mode, listening on port 8000.  
This can be configured by modifying the `cromwellUrl` parameter in `application.conf`.

You can get a build of your current cromwell code with [these instructions](Building.md).
The server can be run with `java -jar <Cromwell JAR> server`, checkout [this page](../CommandLine.md) 
for more detailed instructions. 
You can now run the tests from another terminal.

## Running

There are two ways to invoke the integration tests:

* `sbt "centaur / IntegrationTest / test"` - compiles Centaur and runs all tests via sbt directly. Tests are expected to be in the `centaur/src/main/standardTestCases` directory. This can be changed by modifying `reference.conf`.

* `src/ci/bin/testCentaurLocal.sh` - runs the same tests using the continuous integration pipeline configuration.

* Tests that require different Cromwell and Centaur configurations can be invoked by calling the various scripts in `src/ci/bin`.

### Tags

All tests are tagged with their name and their TESTFORMAT, and also any custom tags specified in the `.test` file.

Tag names are all lower case, so a test named "tagFoo" has a tag "tagfoo".

To run only those tests which have been tagged with a specified tag `tagFoo`:
```
sbt "centaur / IntegrationTest / testOnly * -- -n tagfoo"
```

Or to instead exclude all tests which have been tagged with a specified tag `tagFoo`:
```
sbt "centaur / IntegrationTest / testOnly * -- -l tagfoo"
```

## Adding custom tests

You can add your own tests to the test suite by adding `-Dcentaur.optionalTestPath=DIR` on your sbt invocation, 
e.g. `sbt -Dcentaur.optionalTestPath=/some/path/to/tests test`. The value of `DIR` is expected to be a directory
which contains one or more test case files.
 
The same result can be achieved more permanently by adding the custom directory into the `application.conf` file directly: 
```
centaur {
  optionalTestPath = "/some/path/to/tests"
}
```

## Defining test cases

Each test case file is a HOCON file with the following structure:
```
name: NAME  // Required: Name of the test
testFormat: TESTFORMAT // Required: One of WorkflowSuccessTest, WorkflowFailureTest, runtwiceexpectingcallcaching
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
  // Examples:
  // failures is a list, the first entry (0) might be the error you are looking for. If multiple errors are expected the entire list can be checked. 
  // It has a "message" and a "causedBy" field.
  "failures.0.message": "Cromwell senses you did not use WomTool validate."
  "failures.0.causedBy": "BetweenKeyboardAndChairException"
}

filesystemcheck: "local" // possible values: "local", "gcs". Used in conjunction with outputExpectations to define files we expect to exist after running this workflow.
outputExpectations: {
    "/path/to/my/output/file1": 1
    "/path/to/file/that/should/not/exist": 0
}
```

The tags are optional. If supplied they will allow people to turn on or off this test case by including or excluding tags when running (see above).

The `basePath` field is optional, but if supplied all paths will be resolved from that directory. If it is not supplied, all paths will be resolved from the directory the test case file is in.

The `testFormat` field can be one of the following, case insensitive:
* `workflowsuccess`: The workflow being supplied is expected to successfully complete
* `workflowfailure`: The workflow being supplied is expected to fail

The `metadata` is optional. If supplied, Centaur will retrieve the metadata from the successfully completed workflow and compare the values retrieved to those supplied. At the moment the only fields supported are strings, numbers and booleans.

You can find which metadata is recorded by running a workflow ```java -jar <Cromwell JAR> run -m metadata.json my_workflow.wdl```.
This will save the metadata in `metadata.json`.

For any metadata values or outputExpectations which require workflow ID (i.e, file paths), use `<<UUID>>` as a placeholder instead. For example:
* `"calls.hello.hello.stdout": "gs://google-project/jes/root/wdl/<<UUID>>/call-task/task-stdout.log"`

In case the absolute path the cromwell root is used (for example: `/home/my_user/projects/cromwell/cromwell-executions`)
 you can use `<<WORKFLOW_ROOT>>` as a replacement. 
* `"calls.hello.hello.exit_code": "<<WORKFLOW_ROOT>>/call-hello/execution/exit_code"`

In case testing of the caching is required `<<CACHE_HIT_UUID>>` can be used. 
The testFormat should be `runtwiceexpectingcallcaching`.


## Centaur Test Types
Both Cromwell and Centaur require configuration files in order to correctly build and test various parts of the system. Because of this, we divide
our tests into groups depending on which configuration files they require. Below is the current matrix of configuration files and test source directories.

## Upgrade / Horicromtal / etc.

| CI Test Type                  | Cromwell Config                                                  | Centaur Config                                         |
|-------------------------------|------------------------------------------------------------------|--------------------------------------------------------|
| Engine Upgrade                | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| Horicromtal                   | `papi_[v2beta or v2alpha1]_horicromtal_application.conf`**       | `centaur_application_`<br>`horicromtal.conf`           |
| Horicromtal<br>Engine Upgrade | `papi_v2beta_application.conf`**                                 | `centaur_application_`<br>`horicromtal_no_assert.conf` |
| PAPI Upgrade                  | `papi_v1_v2alpha1_upgrade_application.conf`**                    | `centaur_application.conf`*                            |
| Papi Upgrade<br>New Workflows | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| Azure Blob                    | `centaur_blob_test.conf`**                                       | `centaur_application.conf`*                            |
| WDL Upgrade                   | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |
| (other)                       | `(backend)_application.conf`                                     | `centaur_application.conf`*                            |

| CI Test Type                  | ScalaTest Spec              | Test Directory                      |
|-------------------------------|-----------------------------|-------------------------------------|
| Engine Upgrade                | `EngineUpgradeTestCaseSpec` | `engineUpgradeTestCases`            |
| Horicromtal                   | `CentaurTestSuite`          | `standardTestCases`***              |
| Horicromtal<br>Engine Upgrade | `EngineUpgradeTestCaseSpec` | `engineUpgradeTestCases`***         |
| PAPI Upgrade                  | `PapiUpgradeTestCaseSpec`   | `papiUpgradeTestCases`              |
| PAPI Upgrade<br>New Workflows | `CentaurTestSuite`          | `papiUpgradeNewWorkflowsTestCases`  |
| Azure Blob                    | `CentaurTestSuite      `    | `azureBlobTestCases`                |
| (other)                       | `CentaurTestSuite`          | `standardTestCases`                 |

<small>
\* Centaur Config always uses `centaur_application.conf` except when overridden with `papi_v2alpha1_centaur_application.conf`
or `papi_v2beta_centaur_application.conf`
  ([48 preview link](https://github.com/broadinstitute/cromwell/blob/a7d0601/src/ci/bin/test.inc.sh#L455-L457))  
\*\* Cromwell Config overrides
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/src/ci/bin/test.inc.sh#L213-L221))  
\*\*\* Test Directory overrides
  ([47 link](https://github.com/broadinstitute/cromwell/blob/47/src/ci/bin/test.inc.sh#L440-L449))
</small>

- Engine Upgrade: Retrieves the [Cromwell Version](https://github.com/broadinstitute/cromwell/blob/47/project/Version.scala#L8) then retrieves the previous jar/docker-image from DockerHub. Centaur starts with the prior version, then restarts with the compiled source code.
- Horicromtal: Runs a [docker-compose](https://github.com/broadinstitute/cromwell/blob/47/src/ci/docker-compose/docker-compose-horicromtal.yml) with:
    1. db-mstr: started first
    2. sum-back: runs summarizer
    3. front-back: exposes HTTP
- Horicromtal Engine Upgrade: Combination of Horicromtal and Engine Upgrade
- PAPI Upgrade: Tests run with an older version of Papi and upon restart use a newer version of Papi
- PAPI Upgrade New Workflows: Test definition [does not run any tests](https://travis-ci.org/broadinstitute/cromwell/jobs/475378412)
- WDL Upgrade: Upgrades WDL from draft-2 to 1.0 before testing
- (other): Runs `*.test` files listing the configured backend names

## RDBMS

| Backend | MySQL  | PostgreSQL  | MariaDB  |
|---------|:------:|:-----------:|:--------:|
| AWS     |   ✅   |             |          |
| Local   |   ✅   |      ✅     |          |
| PAPI V2 |   ✅   |             |    ⭕    |
| SLURM   |   ✅   |             |          |
| TES     |   ✅   |             |          |

<small>
⭕ Tests Horicromtal Engine Upgrade versus standard Centaur suite
</small>

All backends run against MySQL. The Local backend also test PostgreSQL, allowing contributors ensure WDLs work with PostgreSQL. MariaDB is tested on a specialized upgrade, where the MySQL connector client is used first, and the MariaDB client is used after restart.
