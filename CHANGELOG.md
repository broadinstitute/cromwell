# Cromwell Change Log

## 41 Release Notes

### Workflow Options

* It is now possible to supply custom `google-labels` in [workflow options](https://cromwell.readthedocs.io/en/stable/wf_options/Google/).

## 40 Release Notes

### Config Changes

#### Cromwell ID in instrumentation path

When set, the configuration value of `system.cromwell_id` will be prepended to StatsD metrics. More info [here](https://cromwell.readthedocs.io/en/stable/developers/Instrumentation/).

#### HealthMonitor Configuration

The HealthMonitor configuration has been refactored to provide a simpler interface:
* You no longer need to specify a monitor class in your `cromwell.conf` as this will now be inherited from the `reference.conf` value.
* You can now opt-in and opt-out of any combination of status monitors.
* The PAPI backends to monitor can now be listed in a single field.

##### Upgrading

You are no longer tied to the previous preset combinations of health checks. However if you just want to carry forward
the exact same set of health checks, you can use one of the following standard recipes:

###### From default, or `NoopHealthMonitorActor`:
If you're currently using the (default) NoopHealthMonitorActor, no action is required.

###### From `StandardHealthMonitorServiceActor`:
If you're currently using the `StandardHealthMonitorServiceActor`, replace this stanza:
```
services {
    HealthMonitor {
        class = "cromwell.services.healthmonitor.impl.standard.StandardHealthMonitorServiceActor"
    }
}
``` 
With this one:
```
services {
    HealthMonitor {
        config {
            check-dockerhub: true
            check-engine-database: true
        }
    }
}
``` 
###### From `WorkbenchHealthMonitorServiceActor`:
Replace this stanza:
```
services {
    HealthMonitor {
        class = "cromwell.services.healthmonitor.impl.workbench.WorkbenchHealthMonitorServiceActor"

        config {
            papi-backend-name = PAPIv1
            papi-v2-backend-name = PAPIv2

            google-auth-name = service-account
            gcs-bucket-to-check = "cromwell-ping-me-bucket"
        }
    }
}
``` 
With this one:
```
services {
    HealthMonitor {
        config {
            check-dockerhub: true
            check-engine-database: true
            check-gcs: true
            check-papi-backends: [PAPIv1, PAPIv2]

            google-auth-name = service-account
            gcs-bucket-to-check = "cromwell-ping-me-bucket"
    }
  }
}
``` 


### Bug fixes

#### WDL 1.0 strings can contain escaped quotes

For example, the statement `String s = "\""` is now supported, whereas previously it produced a syntax error.

#### Empty call blocks in WDL 1.0

Cromwell's WDL 1.0 implementation now allows empty call blocks, e.g. `call task_with_no_inputs {}`. This brings 1.0 in line with draft-2, which has always supported this syntax.

#### Packed CWL bugfix

Fixed a bug that caused an error like `Custom type was referred to but not found` to be issued when using an imported type as a `SchemaDefRequirement` in packed CWL.

## 39 Release Notes

### Cromwell ID changes

When set, the configuration value of `system.cromwell_id` will now have a random suffix appended, unless the
configuration key `system.cromwell_id_random_suffix` is set to `false`.

The generated id also appears more places in the logs, including when picking up workflows from the database and during
shutdown.

### Bug fixes

#### Format fix for `write_map()` 

Fixed an issue that caused the `write_map()` function in Cromwell's WDL 1.0 implementation to produce output in the wrong format. Specifically, the output's rows and columns were swapped. WDL draft-2 was not affected.
  
Incorrect `write_map()` output in Cromwell 38 and earlier:
```
key1    key2    key3
value1  value2  value3
```
Corrected `write_map()` output in Cromwell 39 and later:
```
key1  value1
key2  value2
key3  value3
```

## 38 Release Notes

### HPC paths with Docker

The `ConfigBackendLifecycleActorFactory` path variables `script`, `out` and `err` are now consistent when running with
and without docker. Similarly, when killing a docker task the `kill-docker` configuration key is now used instead of
`kill`. For more information see the [online documentation](https://cromwell.readthedocs.io/en/stable/backends/SGE/).

### No-op Health Monitor is now the default health monitor

Previous versions of Cromwell defaulted to using a health monitor service that checked Docker Hub and engine database status.
Neither check was useful if the `status` endpoint was never consulted as is likely the case in most deployments. Cromwell 38
now defaults to a `NoopHealthMonitorServiceActor` which does nothing. The previous health service implementation is still
available as `StandardHealthMonitorServiceActor`.

### Bug fixes
- Fixed an issue that could cause Cromwell to consume disk space unnecessarily when using zipped dependencies

#### HTTP responses

- When returning errors as json the `Content-Type` header is set to `application/json`.

## 37 Release Notes

### Docker

- Adds support for retrieving docker digests of asia.gcr.io images
- Adds configuration settings for docker digest lookups. See the `docker` section of the `reference.conf` for more information 
- Attempt to automatically adjust the boot disk size on the Google Cloud Backend (version 2) if the size of the image is greater than the default disk size or the required disk size in the runtime attributes.
Only works for registries that support the version 2 of the manifest schema (https://docs.docker.com/registry/spec/manifest-v2-2/)
At this date (12/09/18) this includes GCR and Dockerhub.

### Added new call cache path+modtime hashing strategy.

Call caching hashes with this new strategy are based on the path and the last modified time of the file.

### Instance independent abort

For multi-instance Cromwell deployments sharing a single database, earlier versions of Cromwell required abort
requests to be sent specifically to the instance that was running the targeted workflow. Cromwell 37 now
allows abort commands to be sent to any Cromwell instance in the shared-database deployment. Configuration details
[here](https://cromwell.readthedocs.io/en/develop/Configuring/abort-configuration).


### Call cache blacklisting

The Google Pipelines API (PAPI) version 1 and 2 backends now offer the option of call cache blacklisting on a per-bucket basis.
More info [here](http://cromwell.readthedocs.io/en/develop/CallCaching/#call-cache-copy-authorization-failure-prefix-blacklisting).

### WDL

- All memory units in WDL are now treated as base-2.
For instance `1 KB == 1 KiB == 1024 Bytes`.

### Backend name for call caching purposes

Previous versions of Cromwell incorporated the name of the backend on which a call was run into the call cache hashes generated for that call.
Unfortunately this made it impossible to change the name of a backend without losing all previously run calls as potential cache hits.
Cromwell 37 introduces the `name-for-call-caching-purposes` backend configuration option as a means of decoupling the backend name from the
value used for the backend name for call caching purposes.

### CWL

Support `InputResourceRequirement` hint

### Changing configuration options

#### Logging Token Distribution

In cases where its not obvious why jobs are queued in Cromwell, you can enable logging for the Job Execution Token Dispenser, using
the `system.hog-safety.token-log-interval-seconds` configuration value.

The default, `0`, means that no logging will occur. 

#### HTTP Filesystem

- The HTTP filesystem is now enabled for engine use by default. To continue without an HTTP filesystem, you can add the 
following content into the appropriate stanza of your configuration file:
```
engine {
  filesystems {
    http { 
      enabled: false 
    }
  }
}
``` 
- When the value `exit-code-timeout-seconds` is set, `check-alive` command is now only called once every timeout interval instead of each poll.

### Beta preview of new Womtool `/describe` endpoint

This new endpoint brings the functionality of Womtool to the world of web services. Submit workflows for validation and receive a JSON description in response.

The endpoint is still undergoing heavy development and should not be used in production. The final version will ship in a future release of Cromwell; watch this space.   

### Bug fixes

- Fixed a regression in Cromwell 36 that could cause operations on empty arrays to fail with a spurious type error (closes [#4318](https://github.com/broadinstitute/cromwell/issues/4318))

#### Abort On Hold Workflows

On Hold workflows may now be aborted.

#### Command fixes for AWS and TES

The AWS and TES backends can now handle calls that generate longer command lines. Like the other
backends, commands scripts are first written to a file, the file is downloaded to the execution
host, and then the localized script is run.

Also fixed are AWS `command {}` blocks that use `|` at the start of a line. For example:

```
command {
  echo hello world \
  | cat
}
```

## 36 Release Notes

### Extra configuration options

The value `exit-code-timeout-seconds` can now set in a backend configuration.
Details [here](https://cromwell.readthedocs.io/en/develop/backends/HPC/#exit-code-timeout)

### [AWS S3 file transfers are now encrypted](https://github.com/broadinstitute/cromwell/pull/4264)

### Bug fixes

#### Metadata Request Coalescing

Coalesce metadata requests to eliminate expensive and redundant queries and metadata construction.

#### Eliminate redundant SFS logging and metadata 

Eliminate superfluous logging and metadata publishing in the shared filesystem backend on poll intervals where there was not a state change.

#### AWS region configuration respected throughout

Previously US-EAST-1 was hardcoded in places.

## 35 Release Notes

### Submit workflow using URL

Cromwell now allows for a user to submit the URL pointing to workflow file to run a workflow.
More details on how to use it in: 
- `Server` mode can be found [here](https://cromwell.readthedocs.io/en/develop/api/RESTAPI/).
- `Run` mode can be found [here](https://cromwell.readthedocs.io/en/develop/CommandLine/#run).

### Languages

- Added an opt-in namespace cache for the WDL Draft 2 language factory. Please see the Cromwell example configuration for details. NOTE: if upgrading from a hotfix version of Cromwell
that relied upon this cache, the cache is now opt-in and must be turned on explicitly in config.
- To maintain conformance with the OpenWDL spec, Cromwell drops support for the `version draft-3` identifier in this release. In the rare case where end users may have been using `version draft-3`, `version 1.0` is a drop-in replacement with no effect on functionality.

### HTTP Workflow Inputs for Shared File System and Google Pipelines API Version 2 Backends

`http` and `https` workflow inputs are now supported for shared filesystem and Google Pipelines API (PAPI) version 2
backends. Configuration details are described [here](http://cromwell.readthedocs.io/en/develop/filesystems/HTTP).

### Call cache hint support

More efficient cache hit copying in multi-user environments is now supported through the `call_cache_hit_path_prefixes` workflow option.
Details [here](http://cromwell.readthedocs.io/en/develop/CallCaching/#call-cache-hit-path-prefixes)

### Root workflow level file hash caching support

Cromwell now offers the ability to cache file hashes on a root workflow level basis, details [here](http://cromwell.readthedocs.io/en/develop/CallCaching/#file-hash-caching).

### Extra configuration options

The value `dockerRoot` can now be set in a backend configuration. 
This will set the execution folder in the container (default: `/cromwell-executions`).

### Bug Fixes

#### API
- The `releaseHold` endpoint will now return `404 Not Found` for an unrecognized workflow ID and `400 Bad Request` for a malformed or invalid workflow ID.

#### Languages

- Fixed a bug that allowed values to be "auto-boxed" into a single-element `Array` of that type, which is not allowed in the WDL spec (Closes [#3478](https://github.com/broadinstitute/cromwell/issues/3478)).

#### PAPI version 1

- Restored standard output and error streaming for jobs.

## 34 Release Notes

### Query API

* Fixes a bug which stopped `includeSubworkflow=false` from paging correctly and subworkflows from being discounted correctly from `totalResultsCount`.
* Query results will now be returned in reverse chronological order, with the most-recently submitted workflows returned first.

### Requester Pays on GCS

Access of Google Cloud Storage buckets with Requester Pays enabled is now supported.
Please read the [relevant documentation](http://cromwell.readthedocs.io/en/develop/filesystems/GoogleCloudStorage#requester-pays) for information on how to enable it and the consequences.

### Private Docker Support on Pipelines API v2

Support for private Docker Hub images is now included in the Google Pipelines API v2 backend. PAPI v2 private Docker support is
equivalent to that in PAPI v1 but the configuration differs, please see
[Docker configuration](http://cromwell.readthedocs.io/en/develop/filesystems/Google#Docker) for more details.

### Updated MySQL client with 8.0 support

Updated the MySQL connector client from `5.1.42` to `5.1.46` which adds support for connecting to MySQL 8.0. See the
documentation on [Changes in MySQL Connector/J](https://dev.mysql.com/doc/relnotes/connector-j/5.1/en/news-5-1.html) for
more information.

## 33 Release Notes

### Query endpoint

#### Exclude workflows based on Labels

This gives the ability to **filter out** workflows based on labels. Two new parameters called `excludeLabelAnd` and `excludeLabelOr` can be used for this purpose.
More details on how to use them can be found [here](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/).

#### Include/Exclude subworkflows

Cromwell now supports excluding subworkflows from workflow query results using the `includeSubworkflows` parameter. By default they are included in the results.
More information can be found at [REST API](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/).

#### Query workflows by Submission time

Cromwell now supports querying workflows by submission time. This will help find workflows that are submitted but not started yet (i.e. workflows which are
in On Hold state). More information can be found [here](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/).

#### Submission time in Workflow Query Response

Submission time of a workflow is now included in WorkflowQueryResult, which is part of the response for workflow query.

### File Localization (NIO) Hint

Cromwell now allows tasks in WDL 1.0 can now specify an optimization in their `parameter_meta` that some `File` inputs do not need to be localized for the task to run successfully.
Full details are available in the [documentation page for this optimization](http://cromwell.readthedocs.io/en/develop/optimizations/FileLocalization).

### Bug Fixes

Workflows which are in 'On Hold' state can now be fetched using the query endpoint.

## 32 Release Notes

### Backends

#### Pipelines API V2
Initial support for Google [Pipelines API version 2](https://cloud.google.com/genomics/reference/rest/).
Expect feature parity except for private dockerhub images which are not supported at the moment, but will be in the near future.
Additionally, the "refresh token" authentication mode is **NOT** supported on PAPI V2.

In addition, the following changes are to be expected:
* Error messages for failed jobs might differ from V1
* The Pipelines API log file content might differ from V1

**Important (If you're running Cromwell with a Google backend, read this)**:
The `actor-factory` value for the google backend (`cromwell.backend.impl.jes.JesBackendLifecycleActorFactory`) is being deprecated.
Please update your configuration accordingly.

| PAPI Version |                                 actor-factory                                |
|--------------|:----------------------------------------------------------------------------:|
|      V1      | cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory |
|      V2      | cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory |

If you don't update the `actor-factory` value, you'll get a deprecation warning in the logs, and Cromwell will default back to **PAPI V1**

### Task Retries
Cromwell now supports retrying failed tasks up to a specified count by declaring a value for the [maxRetries](RuntimeAttributes.md#maxRetries) key through the WDL runtime attributes.

### Labels
* Cromwell has removed most of the formatting restrictions from custom labels. Please check the [README](README.md#label-format) for more detailed documentation.
* Custom labels won't be submitted to Google backend as they are now decoupled from Google's default labels.
* Cromwell now publishes the labels as soon as the workflow is submitted (whether started or on hold). If the labels are invalid, the workflow will not be submitted and request will fail.

### Scala 2.11 Removed
From version 32 onwards we will no longer be publishing build artifacts compatible with Scala 2.11. 

* If you don't import the classes into your own scala project then this should have no impact on you.
* If you **are** importing the classes into your own scala project, make sure you are using Scala 2.12.

### Input Validation
Cromwell can now validate that your inputs files do not supply inputs with no impact on the workflow. Strict validation will be disabled by default in WDL draft 2 and CWL but enabled in WDL draft 3. See the 'Language Factory Config' below for details.

### Language Factory Config
All language factories can now be configured on a per-language-version basis. All languages and versions will support the following options:
* `enabled`: Defaults to `true`. Set to `false` to disallow workflows of this language and version.
* `strict-validation`: Defaults to `true` for WDL draft 3 and `false` for WDL draft 2 and CWL. Specifies whether workflows fail if the inputs JSON (or YAML) file contains values which the workflow did not ask for (and will therefore have no effect). Additional strict checks may be added in the future.

### API

* More accurately returns 503 instead of 500 when Cromwell can not respond in a timely manner
* Cromwell now allows a user to submit a workflow but in a state where it will not automatically be picked up for execution. This new state is called 'On Hold'. To do this you need to set the parameter workflowOnHold to true while submitting the workflow.
* API end point 'releaseHold' will allow the user to send a signal to Cromwell to allow a workflow to be startable, at which point it will be picked up by normal execution schemes.

### GPU

The PAPI backend now supports specifying GPU through WDL runtime attributes:

```wdl
runtime {
    gpuType: "nvidia-tesla-k80"
    gpuCount: 2
    zones: ["us-central1-c"]
}
```

The two types of GPU supported are `nvidia-tesla-k80` and `nvidia-tesla-p100`

**Important**: Before adding a GPU, make sure it is available in the zone the job is running in: https://cloud.google.com/compute/docs/gpus/

### Job Shell

Cromwell now allows for system-wide or per-backend job shell configuration for running user commands rather than always
using the default `/bin/bash`. To set the job shell on a system-wide basis use the configuration key `system.job-shell` or on a
per-backend basis with `<config-key-for-backend>.job-shell`. For example:

```
# system-wide setting, all backends get this
-Dsystem.job-shell=/bin/sh
```

```
# override for just the Local backend
-Dbackend.providers.Local.config.job-shell=/bin/sh
```

For the Config backend the value of the job shell will be available in the `${job_shell}` variable. See Cromwell's `reference.conf` for an example
of how this is used for the default configuration of the `Local` backend.

### Bug Fixes

The imports zip no longer unpacks a single (arbitrary) internal directory if it finds one (or more). Instead, import statements should now be made relative to the base of the import zip root.

#### Reverting Custom Labels

Reverting to a prior custom label value now works.

["Retrieves the current labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#retrieves-the-current-labels-for-a-workflow)
will return the most recently summarized custom label value.

The above endpoint may still return the prior value for a short period of time after using
["Updated labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#update-labels-for-a-workflow)
until the background metadata summary process completes.

#### Deleting Duplicate Custom Label Rows

If you never used the REST API to revert a custom label back to a prior value you will not be affected. This only applies to workflows previously updated using
["Updated labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#update-labels-for-a-workflow).

The database table storing custom labels will delete duplicate rows for any workflow label key. For efficiency purposes
the values are not regenerated automatically from the potentially large metadata table.

In rare cases where one tried to revert to a prior custom label value you may continue to see different results
depending on the REST API used. After the database update
["Retrieves the current labels for a workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#retrieves-the-current-labels-for-a-workflow)
will return the most-recent-unique value while
["Get workflow and call-level metadata for a specified workflow"](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#get-workflow-and-call-level-metadata-for-a-specified-workflow)
will return the up-to-date value. For example, if one previously updated a value from `"value-1"` > `"value-2"` >
`"value-3"` > `"value-2"` then the former REST API will return `value-3` while the latter will return `value-2`.

#### Workflow options `google_project` output in metadata

Workflow metadata for jobs run on a Google Pipelines API backend will report the `google_project` specified via a
[workflow options json](http://cromwell.readthedocs.io/en/develop/wf_options/Google/#google-pipelines-api-workflow-options).

## 31 Release Notes

* **Cromwell server**  
The Cromwell server source code is now located under `server/src`. `sbt assembly` will build the runnable Cromwell JAR in 
`server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`.

* **Robustness**
    + The rate at which jobs are being started can now be controlled using the `system.job-rate-control` configuration stanza.  
    + A load controller service has been added to allow Cromwell to self-monitor and adjust its load accordingly.
The load controller is currently a simple on/off switch controlling the job start rate. It gathers metrics from different parts of the system
to inform its decision to stop the creation of jobs.
You can find relevant configuration in the `services.LoadController` section of the `cromwell.examples.conf` file,
as well as in the `load-control` section in `reference.conf`.
The load level of the monitored sub-systems are instrumented and can be found under the `cromwell.load` statsD path.
    + The statsD metrics have been re-shuffled a bit. If you had a dashboard you might find that you need to update it.
Changes include: 
        + Removed artificially inserted "count" and "timing" the path
        + Added a `load` section
        + Metrics were prefixed twice with `cromwell` (`cromwell.cromwell.my_metric`), now they're only prefixed once
        + Added `processed` and `queue` metrics under various metrics monitoring the throughput and amount of queued work respectively
        + Added a memory metric representing an estimation of the free memory Cromwell thinks it has left

* Added a configuration option under `docker.hash-lookup.enabled` to disable docker hash lookup.
 Disabling it will also disable call caching for jobs with floating docker tags.
 
* **API**    
    + Updated the `/query` response to include the total number of query results returned. See [here](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#workflowqueryresponse) for more information.

## 30.1 Release Notes

* A set of bug fixes following the migration of Cromwell to WOM (the Workflow Object Model) in version 30.

## 30 Release Notes

### Breaking changes

* The `customLabels` form field for workflow submission has been renamed to `labels`.

### Other changes

* **New Cromwell documentation**  
Our documentation has moved from our [README](https://github.com/broadinstitute/cromwell/blob/29_hotfix/README.md) to a new website: [Cromwell Documentation](http://cromwell.readthedocs.io/en/develop/). There are new [Tutorials](http://cromwell.readthedocs.io/en/develop/tutorials/FiveMinuteIntro/) and much of the documentation has been re-written. The source files are in the [/docs](https://github.com/broadinstitute/cromwell/tree/develop/docs) directory.

* **API**  
    + Cromwell now supports input files in the yaml format (JSON format is still supported).
    + Added a [GET version for the `labels` endpoint](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#retrieves-the-current-labels-for-a-workflow) which will return current labels for a workflow.

* **Database**  
You have the option of storing the metadata in a separate SQL database than the database containing the internal engine
data. When switching connection information for an existing database containing historical data, the tables
should be manually replicated from one database instance to another using the tools appropriate for your specific
database types. Cromwell will not move any existing data automatically. This feature should be considered experimental
and likely to change in the future. See the [Database Documentation](https://cromwell.readthedocs.io/en/develop/Configuring/#database) or the `database` section in
[cromwell.examples.conf](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.examples.conf) for more
information.

* **StatsD**  
Added initial support for StatsD instrumentation. See the [Instrumentation Documentation](https://cromwell.readthedocs.io/en/develop/Instrumentation) for details on how to use it.

* **User Service Account auth mode for Google**  
Added a new authentication mode for [Google Cloud Platform](https://cromwell.readthedocs.io/en/develop/backends/Google) which will allow a user to supply the JSON key file in their workflow options to allow for per-workflow authentication via service account. This is analogous to the previously existing refresh token authentication scheme. As with the refresh token scheme it is encouraged that the **user_service_account_json** workflow option field is added to the **encrypted-fields** list in the configuration.

* **Bugfixes**  
Abort of Dockerized tasks on the Local backend should now work as expected. Cromwell uses `docker kill` to kill the Docker container.

## 29 Release Notes

### Breaking Changes

* **Command line**  
In preparation for supporting CWL scripts (yes, you read that right!), we have extensively revised the Command Line in Cromwell 29. For more details about the usage changes please see the [README](https://github.com/broadinstitute/cromwell#command-line-usage). And stay tuned to the [WDL/Cromwell blog](https://software.broadinstitute.org/wdl/blog) over the next couple of months for more news about CWL.

* **Request timeouts**   
Cromwell now returns more specific `503 Service Unavailable` error codes on request timeouts, rather than the more generic `500 Internal Server Error`. The response for a request timeout will now be plain text, rather than a JSON format.

* **Metadata endpoint**  
The response from the metadata endpoint can be quite large depending on your workflow. You can now opt-in to have Cromwell gzip your metadata file, in order to reduce file size, by sending the `Accept-Encoding: gzip` header. The default behavior now does not gzip encode responses.

* **Engine endpoints**  
Previously the engine endpoints were available under `/api/engine` but now the endpoints are under `/engine` so they don't require authentication. Workflow endpoints are still available under `/api/workflows`. We also deprecated the setting `api.routeUnwrapped` as a part of this internal consistency effort.

* **Call caching diff**  
We updated the response format of the [callcaching/diff](https://github.com/broadinstitute/cromwell#get-apiworkflowsversioncallcachingdiff) endpoint.

### Other changes

* **Cromwell server**  
When running in server mode, Cromwell now attempts to gracefully shutdown after receiving a `SIGINT` (`Ctrl-C`) or `SIGTERM` (`kill`) signal. This means that Cromwell waits for all pending database writes before exiting, as long as you include `application.conf` at the top of your config file. You can find detailed information about how to configure this feature in the [Cromwell Wiki](https://github.com/broadinstitute/cromwell/wiki/DevZone#graceful-server-shutdown).

* **Concurrent jobs**  
You can now limit the number of concurrent jobs for any backend. Previously this was only possible in some backend implementations. Please see the [README](https://github.com/broadinstitute/cromwell#backend-job-limits) for details.

### WDL

* **Optional WDL variables**  
Empty optional WDL values are now rendered as the `null` JSON value instead of the JSON string `"null"` in the metadata and output endpoints. You do not need to migrate previous workflows. Workflows run on Cromwell 28 and prior will still render empty values as `"null"`.

* **Empty WDL variables**  
Cromwell now accepts `null` JSON values in the input file and coerces them as an empty WDL value. WDL variables must be declared optional in order to be supplied with a `null` JSON value.

input.json
```json
{
    "null_input_values.maybeString": null,
    "null_input_values.arrayOfMaybeInts": [1, 2, null, 4]
}
```

workflow.wdl
```
workflow null_input_values {
    String? maybeString
    Array[Int?] arrayOfMaybeInts
}
```

## 28

### Bug Fixes

#### WDL write_* functions add a final newline

The following WDL functions now add a newline after the final line of output (the previous behavior of not adding this
newline was inadvertent):
- `write_lines`
- `write_map`
- `write_object`
- `write_objects`
- `write_tsv`

For example:

```
task writer {
  Array[String] a = ["foo", "bar"]
  command {
    # used to output: "foo\nbar"
    # now outputs: "foo\nbar\n"
    cat write_lines(a)
  }
}
```

#### `ContinueWhilePossible`

A workflow utilizing the WorkflowFailureMode Workflow Option `ContinueWhilePossible` will now successfully reach a terminal state once all runnable jobs have completed.
#### `FailOnStderr` 
When `FailOnStderr` is set to false, Cromwell no longer checks for the existence of a stderr file for that task. 

### WDL Functions

#### New functions: floor, ceil and round:

Enables the `floor`, `ceil` and `round` functions in WDL to convert floating point numbers to integers.

For example we can now use the size of an input file to influence the amount of memory the task is given. In the example below a 500MB input file will result in a request for a VM with 2GB of memory:

```
task foo {
    File in_file
    command { ... }
    runtime {
      docker: "..."
      memory: ceil(size(in_file)) * 4 
    }
}
```

### Call Caching

* Hash values calculated by Cromwell for a call when call caching is enabled are now published to the metadata.
It is published even if the call failed. However if the call is attempted multiple times (because it has been preempted for example),
since hash values are strictly identical for all attempts, they will only be published in the last attempt section of the metadata for this call.
If the hashes fail to be calculated, the reason is indicated in a `hashFailures` field in the `callCaching` section of the call metadata.
*Important*: Hashes are not retroactively published to the metadata. Which means only workflows run on Cromwell 28+ will have hashes in their metadata.

See the [README](https://github.com/broadinstitute/cromwell#get-apiworkflowsversionidmetadata) for an example metadata response.

* New endpoint returning the hash differential for 2 calls. 

`GET /api/workflows/:version/callcaching/diff`

See the [README](https://github.com/broadinstitute/cromwell#get-apiworkflowsversioncallcachingdiff) for more details.

### Workflow Submission

* The workflow submission parameters `wdlSource` and `wdlDependencies` have been deprecated in favor of `workflowSource` and
`workflowDependencies` respectively.  The older names are still supported in Cromwell 28 with deprecation warnings but will
be removed in a future version of Cromwell.

### Labels
* A new `/labels` endpoint has been added to update labels for an existing workflow. See the [README](README.md#patch-apiworkflowsversionidlabels) for more information.
* Label formatting requirements have been updated, please check the [README](README.md#label-format) for more detailed documentation.


### JES Backend

The JES backend now supports a `filesystems.gcs.caching.duplication-strategy` configuration entry.
It can be set to specify the desired behavior of Cromwell regarding call outputs when a call finds a hit in the cache.
The default value is `copy` which will copy all output files to the new call directory.
A second value is allowed, `reference`, that will instead point to the original output files, without copying them.


```hocon
filesystems {
  gcs {
    auth = "application-default"
    
    caching {
      duplication-strategy = "reference"
    }
  }
}
```

A placeholder file will be placed in the execution folder of the cached call to explain the absence of output files and point to the location of the original ones.


### Metadata Write Batching

Metadata write batching works the same as in previous versions of Cromwell, but the default batch size has been changed from 1 to 200.  It's possible that 200 is too high in some environments, but 200 is more likely to be an appropriate value
than the previous default.


## 27

### Migration

* Call Caching has been improved in this version of Cromwell, specifically the time needed to determine whether or not a job can be cached
 has drastically decreased. To achieve that the database schema has been modified and a migration is required in order to preserve the pre-existing cached jobs.
 This migration is relatively fast compared to previous migrations. To get an idea of the time needed, look at the size of your `CALL_CACHING_HASH_ENTRY` table.
 As a benchmark, it takes 1 minute for a table with 6 million rows.
 The migration will only be executed on MySQL. Other databases will lose their previous cached jobs.
 In order to run properly on MySQL, **the following flag needs to be adjusted**: https://dev.mysql.com/doc/refman/5.5/en/server-system-variables.html#sysvar_group_concat_max_len
 The following query will give you a minimum to set the group_concat_max_len value to:
 
 ```sql
SELECT MAX(aggregated) as group_concat_max_len FROM
      (
            SELECT cche.CALL_CACHING_ENTRY_ID, SUM(LENGTH(CONCAT(cche.HASH_KEY, cche.HASH_VALUE))) AS aggregated
            FROM CALL_CACHING_HASH_ENTRY cche
            GROUP BY cche.CALL_CACHING_ENTRY_ID
      ) aggregation
 ```

 Here is the SQL command to run to set the group_concat_max_len flag to the proper value:
 
 ```sql
SET GLOBAL group_concat_max_len = value
 ```
 
 Where `value` is replaced with the value you want to set it to.
 
 Note that the migration will fail if the flag is not set properly.
 
### Breaking Changes

* The update to Slick 3.2 requires a database stanza to
[switch](http://slick.lightbend.com/doc/3.2.0/upgrade.html#profiles-vs-drivers) from using `driver` to `profile`.

```hocon
database {
  #driver = "slick.driver.MySQLDriver$" #old
  profile = "slick.jdbc.MySQLProfile$"  #new
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

### Call Caching

Cromwell now supports call caching with floating Docker tags (e.g. `docker: "ubuntu:latest"`). Note it is still considered
a best practice to specify Docker images as hashes where possible, especially for production usages.

Within a single workflow Cromwell will attempt to resolve all floating tags to the same Docker hash, even if Cromwell is restarted
during the execution of a workflow. In call metadata the `docker` runtime attribute is now the same as the
value that actually appeared in the WDL:

```
   "runtimeAttributes": {
     "docker": "ubuntu:latest",
     "failOnStderr": "false",
     "continueOnReturnCode": "0"
   }
```

Previous versions of Cromwell rewrote the `docker` value to the hash of the Docker image.

There is a new call-level metadata value `dockerImageUsed` which captures the hash of the Docker image actually used to
run the call:

```
   "dockerImageUsed": "library/ubuntu@sha256:382452f82a8bbd34443b2c727650af46aced0f94a44463c62a9848133ecb1aa8"
```

### Docker

* The Docker section of the configuration has been slightly reworked 
An option to specify how a Docker hash should be looked up has been added. Two methods are available.
    "local" will try to look for the image on the machine where cromwell is running. If it can't be found, Cromwell will try to `pull` the image and use the hash from the retrieved image.
    "remote" will try to look up the image hash directly on the remote repository where the image is located (Docker Hub and GCR are supported)
Note that the "local" option will require docker to be installed on the machine running cromwell, in order for it to call the docker CLI.
* Adds hash lookup support for public [quay.io](https://quay.io/) images.

### WDL Feature Support
* Added support for the new WDL `basename` function. Allows WDL authors to get just the file name from a File (i.e. removing the directory path)
* Allows coercion of `Map` objects into `Array`s of `Pair`s. This also allows WDL authors to directly scatter over WDL `Map`s.

### Miscellaneous
* Adds support for JSON file format for google service account credentials. As of Cromwell 27, PEM credentials for PAPI are deprecated and support might be removed in a future version.

```
google {

  application-name = "cromwell"

  auths = [
    {
      name = "service-account"
      scheme = "service_account"
      json-file = "/path/to/file.json"
    }
  ]
}
```

### General Changes

* The `/query` endpoint now supports querying by `label`. See the [README](README.md#get-apiworkflowsversionquery) for more information.
* The `read_X` standard library functions limit accepted filesizes.  These differ by type, e.g. read_bool has a smaller limit than read_string.  See reference.conf for default settings.

## 26

### Breaking Changes

* Failure metadata for calls and workflows was being displayed inconsistently, with different formats depending on the originating Cromwell version. Failures will now always present as an array of JSON objects each representing a failure. Each failure will have a message and a causedBy field. The causedBy field will be an array of similar failure objects. An example is given below:

```
failures: [{
  message: "failure1",
  causedBy: [{
    message: "cause1",
    causedBy: []
   }, {
    message: "cause2",
    causedBy: []
  }]
 }, {
  message: "failure2",
  causedBy: []
}]
```

### Additional Upgrade Time

* Upgrading to Cromwell 26 will take additional time due to the migration of failure metadata. Cromwell will automatically run a database query during the upgrade which appears to be roughly linear to the number of rows in the METADATA_ENTRY table. You can estimate upgrade time using the following equation: `time to migrate (in seconds) ~= (rows in METADATA_ENTRY) / 65000` Note that due to differences in hardware and database speed, this is only a rough estimate.

### Config Changes

* Added a configuration option under `system.io` to throttle the number of I/O queries that Cromwell makes, as well as configure retry parameters.
 This is mostly useful for the JES backend and should be updated to match the GCS quota available for the project.
 
```
system.io {
  # Global Throttling - This is mostly useful for GCS and can be adjusted to match
  # the quota availble on the GCS API
  number-of-requests = 100000
  per = 100 seconds
  
  # Number of times an I/O operation should be attempted before giving up and failing it.
  number-of-attempts = 5
}
```

## 25

### External Contributors
* A special thank you to @adamstruck, @antonkulaga and @delocalizer for their contributions to Cromwell.
### Breaking Changes

* Metadata keys for call caching are changed. All call caching keys are now in a `callCaching` stanza. `Call cache read result` has moved here and is now `result`. The `allowResultReuse` and `effectiveCallCachingMode` have moved here. The `hit` boolean is a simple indication of whether or not it was a hit, with no additional information. An example using the new format is:
```
"callCaching": {
  "hit": false,
  "effectiveCallCachingMode": "ReadAndWriteCache",
  "result": "Cache Miss",
  "allowResultReuse": true
}
```

### Config Changes

* Added a field `insert-batch-size` to the `database` stanza which defines how many values from a batch insert will be processed at a time. This value defaults to 2000. 
* Moved the config value `services.MetadataService.metadata-summary-refresh-interval` to `services.MetadataService.config.metadata-summary-refresh-interval`
* Added ability to override the default zone(s) used by JES via the config structure by setting `genomics.default-zones` in the JES configuration
* The cromwell server TCP binding timeout is now configurable via the config key `webservice.binding-timeout`, defaulted
  to the previous value `5s` (five seconds) via the reference.conf.
* For MySQL users, a massive scalability improvement via batched DB writing of internal metadata events. Note that one must add `rewriteBatchedStatements=true` to their JDBC URL in their config in order to take advantage of this

### General Changes

* Cromwell's WDL parser now recognizes empty array literals correctly, e.g. `Array[String] emptyArray = []`.
* Cromwell now applies default labels automatically to JES pipeline runs.
* Added support for new WDL functions:
  * `length: (Array[X]) => Integer` - report the length of the specified array
  * `prefix: (String, Array[X]) => Array[String]` - generate an array consisting of each element of the input array prefixed
     by a specified `String`.  The input array can have elements of any primitive type, the return array will always have
     type `Array[String]`.
  * `defined: (Any) => Boolean` - Will return false if the provided value is an optional that is not defined. Returns true in all other cases.
* Cromwell's Config (Shared Filesystem) backend now supports invocation of commands which run in a Docker image as a non-root user.
  The non-root user could either be the default user for a given Docker image (e.g. specified in a Dockerfile via a `USER` directive),
  or the Config backend could pass an optional `"-u username"` as part of the `submit-docker` command.
* In some cases the SFS backend, used for Local, SGE, etc., coerced `WdlFile` to `WdlString` by using `.toUri`. This
resulted in strings prepended with `file:///path/to/file`. Now absolute file paths will not contain the uri scheme.
* Launch jobs on servers that support the GA4GH Task Execution Schema using the TES backend.
* **Call caching: Cromwell will no longer try to use the cache for WDL tasks that contain a floating docker tag.** 
  Call caching will still behave the same for tasks having a docker image with a specific hash.
  See https://github.com/broadinstitute/cromwell#call-caching-docker-tags for more details. 
* Added docker hash lookup. Cromwell will try to lookup the hash for a docker image with a floating tag, and use that hash when executing the job.
  This will be reflected in the metadata where the docker runtime attribute will contains the hash that was used.
  If Cromwell is unable to lookup the docker hash, the job will be run with the original user defined floating tag.
  Cromwell is currently able to lookup public and private docker hashes for images on Docker Hub and Google Container Engine for job running on the JES backend.
  For other backends, cromwell is able to lookup public docker hashes for Docker Hub and Google Container Engine.
  See https://github.com/broadinstitute/cromwell#call-caching-docker-tags for more details. 

### Database schema changes
* Added CUSTOM_LABELS as a field of WORKFLOW_STORE_ENTRY, to store workflow store entries.

## 24

* When emitting workflow outputs to the Cromwell log only the first 1000 characters per output will be printed
* Added support for conditional (`if`) statements.
* Globs for Shared File System (SFS) backends, such as local or SGE, now use bash globbing instead of Java globbing, consistent with the JES backend.

## 23

* The `meta` and `parameter_meta` blocks are now valid within `workflow` blocks, not just `task`
* The JES backend configuration now has an option `genomics-api-queries-per-100-seconds` to help tune the rate of batch polling against the JES servers. Users with quotas larger than default should make sure to set this value.
* Added an option `call-caching.invalidate-bad-cache-results` (default: `true`). If true, Cromwell will invalidate cached results which have failed to copy as part of a cache hit.
* Timing diagrams and metadata now receive more fine grained workflow states between submission and Running.
* Support for the Pair WDL type (e.g. `Pair[Int, File] floo = (3, "gs://blar/blaz/qlux.txt")`)
* Added support for new WDL functions:
  * `zip: (Array[X], Array[Y]) => Array[Pair[X, Y]]` - align items in the two arrays by index and return them as WDL pairs 
  * `cross: (Array[X], Array[Y]) => Array[Pair[X, Y]]` - create every possible pair from the two input arrays and return them all as WDL pairs
  * `transpose: (Array[Array[X]]) => Array[Array[X]]` compute the matrix transpose for a 2D array. Assumes each inner array has the same length.
* By default, `system.abort-jobs-on-terminate` is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <wdl> <inputs>`.
* Enable WDL imports when running in Single Workflow Runner Mode.
* Both batch and non-batch REST workflow submissions now require a multipart/form-data encoded body.
* Support for sub workflows (see [Annex A](#annex-a---workflow-outputs))
* Enable WDL imports when running in Single Workflow Runner Mode as well as Server Mode
* Support for WDL imports through an additional imports.zip parameter
* Support for sub workflows
* Corrected file globbing in JES to correctly report all generated files. Additionally, file globbing in JES now uses bash-style glob syntax instead of python style glob syntax
* Support declarations as graph nodes
* Added the ability to override the default service account that the compute VM is started with via the configuration option `JES.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`. More details can be found in the README.md
* Fix bugs related to the behavior of Cromwell in Single Workflow Runner Mode. Cromwell will now exit once a workflow completes in Single Workflow Runner Mode. Additionally, when restarting Cromwell in Single Workflow Runner Mode, Cromwell will no longer restart incomplete workflows from a previous session.

### Annex A - Workflow outputs
    
The WDL specification has changed regarding [workflow outputs](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#outputs) to accommodate sub workflows.
This change is backward compatible in terms of runnable WDLs (WDL files using the deprecated workflow outputs syntax will still run the same). 
The only visible change lies in the metadata (as well as the console output in single workflow mode, when workflow outputs are printed out at the end of a successful workflow).

TL;DR Unless you are parsing or manipulating the "key" by which workflow outputs are referenced in the metadata (and/or the console output for single workflow mode), you can skip the following explanation.

*Metadata Response*
```
{
  ...
  outputs {
    "task_output_1": "hello",
    "task_output_2": "world"
            ^
       If you don't manipulate this part of the metadata, then skip this section
  }
}
```

In order to maintain backward compatibility, workflow outputs expressed with the deprecated syntax are "expanded" to the new syntax. Here is an example:

```
task t {
    command {
        #do something
    }
    output {
        String out1 = "hello"
        String out2 = "world"
    }
}
```

```
    workflow old_syntax {
        call t
        output {
            t.*
        }
    }
```

```
    workflow new_syntax {
        call t
        output {
            String wf_out1 = t.out1
            String wf_out2 = t.out2
        }
    }
```

The new syntax allows for type checking of the outputs as well as expressions. It also allows for explicitly naming to the outputs.
The old syntax doesn't give the ability to name workflow outputs. For consistency reasons, Cromwell will generate a "new syntax" workflow output for each task output, and name them.
Their name will be generated using their FQN, which would give 

```
output {
   String w.t.out1 = t.out1
   String w.t.out2 = t.out2
}
```
        
However as the FQN separator is `.`, the name itself cannot contain any `.`. 
For that reason, `.` are replaced with `_` :

*Old syntax expanded to new syntax*
```
output {
   String w_t_out1 = t.out1
   String w_t_out2 = t.out2
}
```

The consequence is that the workflow outputs section of the metadata for `old_syntax` would previously look like 
 
 ```
    outputs {
        "w.t.out1": "hello",
        "w.t.out2": "hello"
    }
 ```
 
but it will now look like 

```
    outputs {
        "w_t_out1": "hello",
        "w_t_out2": "hello"
    }
```

The same applies for the console output of a workflow run in single workflow mode.


## 0.22

* Improved retries for Call Caching and general bug fixes.
* Users will experience better scalability of status polling for Google JES.
* Now there are configurable caching strategies for a SharedFileSystem backend (i.e. Local, SFS) in the backend's stanza:
  See below for detailed descriptions of each configurable key.

```
backend {
  ...
  providers {
    SFS_BackendName {
      actor-factory = ...
      config {
        ...
        filesystems {
          local {
            localization: [
               ...
            ]
            caching {
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]
              # Possible values: file, path
              # "file" will compute an md5 hash of the file content.
              # "path" will compute an md5 hash of the file path. This strategy will only be effective if the duplication-strategy (above) is set to "soft-link",
              # in order to allow for the original file path to be hashed.
              hashing-strategy: "file"

              # When true, will check if a sibling file with the same name and the .md5 extension exists, and if it does, use the content of this file as a hash.
              # If false or the md5 does not exist, will proceed with the above-defined hashing strategy.
              check-sibling-md5: false
            }
```
* Multiple Input JSON files can now be submitted in server mode through the existing submission endpoint: /api/workflows/:version.
    This endpoint accepts a POST request with a multipart/form-data encoded body. You can now include multiple keys for workflow inputs.

        Each key below can contain an optional JSON file of the workflow inputs. A skeleton file can be generated from wdltool using the "inputs" subcommand.
        NOTE: In case of key conflicts between multiple JSON files, higher values of x in workflowInputs_x override lower values. For example, an input
        specified in workflowInputs_3 will override an input with the same name that was given in workflowInputs or workflowInputs_2. Similarly, an input
        specified in workflowInputs_5 will override an input with the same name in any other input file.

        workflowInputs
        workflowInputs_2
        workflowInputs_3
        workflowInputs_4
        workflowInputs_5

* You can now limit the number of concurrent jobs for a backend by specifying the following option in the backend's config stanza:
```
backend {
  ...
  providers {
    BackendName {
      actor-factory = ...
      config {
        concurrent-job-limit = 5
```


## 0.21

* Warning: Significant database updates when you switch from version 0.19 to 0.21 of Cromwell.
  There may be a long wait period for the migration to finish for large databases.
  Please refer to MIGRATION.md for more details.

* There are significant architectural changes related to increases in performance and scaling.

* The biggest user-facing changes from 0.19 to 0.21 are related to the application.conf file, which has been restructured significantly.
The configuration for backends now is all contained within a `backend` stanza, which specifies 1 stanza per name per backend and a default backend, as follows:

```
backend {
    default=Local
    providers {
        Local {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                ... backend specific config ...
            }
        }
        JES {
            actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
            config {
                ... backend specific config ...
            }
        }
        SGE {
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
            config {
                ... backend specific config ...
            }r
        }
    }
}
```
* A new `/stats` endpoint has been added to get workflow and job count for a Cromwell running in server mode.

* Renamed Workflow Options:
   “workflow_log_dir” -> “final_workflow_log_dir”
    “call_logs_dir” -> “final_call_logs_dir”
    “outputs_path” -> “final_workflow_outputs_dir”
    “defaultRuntimeOptions” -> “default_runtime_attributes”

* Timing diagrams endpoint has been updated to include additional state information about jobs.

* Add support for Google Private IPs through `noAddress` runtime attribute. If set to true, the VM will NOT be provided with a public IP address.
*Important*: Your project must be whitelisted in "Google Access for Private IPs Early Access Program". If it's not whitelisted and you set this attribute to true, the task will hang.
  Defaults to `false`.
  e.g:
```
task {
    command {
        echo "I'm private !"
    }

    runtime {
        docker: "ubuntu:latest"
        noAddress: true
    }
}
```

* The Local and the SGE backend have been merged into a generic
Shared File System (SFS) backend. This updated backend can be configured
to work with various other command line dispatchers such as LSF. See the
[README](README.md#sun-gridengine-backend) for more info.

* On the JES and SFS backends, task `command` blocks are now always
passed absolute paths for input `File`s.

* On the SFS backends, the call directory now contains two sub-directories:
    * `inputs` contains all the input files that have been localized for this task (see next below for more details)
    * `execution` contains all other files (script, logs, rc, potential outputs etc...)

* Override the default database configuration by setting the keys
`database.driver`, `database.db.driver`, `database.db.url`, etc.
* Override the default database configuration by setting the keys
`database.driver`, `database.db.driver`, `database.db.url`, etc.

For example:
```
# use a mysql database
database {
  driver = "slick.driver.MySQLDriver$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

## 0.20

* The default per-upload bytes size for GCS is now the minimum 256K
instead of 64M. There is also an undocumented config key
`google.upload-buffer-bytes` that allows adjusting this internal value.

* Updated Docker Hub hash retriever to parse json with [custom media
types](https://github.com/docker/distribution/blob/05b0ab0/docs/spec/manifest-v2-1.md).

* Added a `/batch` submit endpoint that accepts a single wdl with
multiple input files.

* The `/query` endpoint now supports querying by `id`, and submitting
parameters as a HTTP POST.
