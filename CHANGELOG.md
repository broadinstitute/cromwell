# Cromwell Change Log

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
    driver = "com.mysql.jdbc.Driver"
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
    
The WDL specification has changed regarding [workflow outputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs) to accommodate sub workflows.
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
    driver = "com.mysql.jdbc.Driver"
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