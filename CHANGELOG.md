# Cromwell Change Log

## 80 Release Notes

### Direct WES support in Cromwell

Cromwell 80 no longer supports the wes2cromwell project within the Cromwell repository.

In the previous release, 3 Wes2Cromwell endpoints in the Cromwell project were implemented and documented in the Swagger API. Three new endpoints,
located within the wes2cromwell project, will also be moved, implemented, and documented within Cromwell. As a result of this, we can safely remove 
and deprecate the wes2cromwell project from the repo.

Previous endpoints:

| HTTP verb | Endpoint path | Description   |
| --------- | ------------- |---------------|
| GET | /api/ga4gh/wes/v1/service-info | Server info |
| POST | /api/ga4gh/wes/v1/runs/{run_id}/cancel | Abort workflow |
| GET | /api/ga4gh/wes/v1/runs/{run_id}/status | Workflow status |

Newly implemented endpoints:

| HTTP verb | Endpoint path | Description     |
| --------- | ------------- |-----------------|
| GET | /api/ga4gh/wes/v1/runs | List workflows  |
| POST | /api/ga4gh/wes/v1/runs | Submit workflow |
| GET | /api/ga4gh/wes/v1/runs/{run_id} | Workflow details |

### Batch Compute Service backend removed

The BCS backend and OSS filesystem (both of which support Alibaba Cloud) have been removed.

## 79 Release Notes

### Last release with CWL support

Cromwell 79 is the last release with CWL. Support will be removed in Cromwell 80 and above.

CWL will be re-introduced at a later date in the [Terra platform](https://terra.bio/), using a solution other than Cromwell. See the blog post ["Terra’s roadmap to supporting more workflow languages"](https://terra.bio/terras-roadmap-to-supporting-more-workflow-languages/) for details.

| Product                                   | Language | Support               |
|-------------------------------------------|----------|-----------------------|
| Cromwell standalone                       | WDL      | :white_check_mark:    |
| Cromwell standalone                       | CWL      | :x:                    |
| [Terra SaaS platform](https://terra.bio/) | WDL      | :white_check_mark:    |
| [Terra SaaS platform](https://terra.bio/) | CWL      | Future support planned |

### Last release with Alibaba Cloud

The BCS backend and OSS filesystem (both of which support Alibaba Cloud) will be removed in version 80.

### WES endpoints preview

As a means to stay on top of endpoints within our repo, 3 new Workflow Execution Service (WES) endpoints are now documented in the Cromwell Swagger (others to follow as part of later work):

| HTTP verb | Endpoint path | Description   |
| --------- | ------------- |---------------|
| GET | /api/ga4gh/wes/v1/service-info | Server info |
| POST | /api/ga4gh/wes/v1/runs/{run_id}/cancel | Abort workflow |
| GET | /api/ga4gh/wes/v1/runs/{run_id}/status | Workflow status |

### Scala 2.13

Cromwell is now built with Scala version 2.13. This change should not be noticeable to users but may be of interest to developers of Cromwell backend implementations.

### Bug Fixes

 * Fixed a call caching bug in which an invalid cache entry could cause a valid cache entry to be ignored.

## 75 Release Notes

### New `AwaitingCloudQuota` backend status

For Cloud Life Sciences v2beta only.

When a user's GCP project reaches a quota limit, Cromwell continues to submit jobs and Life Sciences acknowledges them as created even if the physical VM cannot yet start. Cromwell now detects this condition in the backend and reports `AwaitingCloudQuota`.

The status is informational and does not require any action. Users wishing to maximize throughput can use `AwaitingCloudQuota` as an indication they should check quota in Cloud Console and request a quota increase from GCP.

`AwaitingCloudQuota` will appear between the `Initializing` and `Running` backend statuses, and will be skipped if not applicable.

Now:

| Status in metadata |Quota normal| Quota delay          | Status meaning                                    |
|--------------------|----|----------------------|---------------------------------------------------|
| `executionStatus`    |`Running`| `Running`            | Job state Cromwell is requesting from the backend |
| `backendStatus`      |`Running`| `AwaitingCloudQuota` | Job state reported by backend                          |

Previously:

| Status in metadata |Quota normal|Quota delay| Status meaning                                            |
|--------------------|----|----|-----------------------------------------------------------|
| `executionStatus`    |`Running`|`Running`| Job state Cromwell is requesting from the backend |
| `backendStatus`      |`Running`|`Running`| Job state reported by backend |

### New 'requestedWorkflowId' API Option

Allows users to choose their own workflow IDs at workflow submission time. 

If supplied for single workflows, this value must be a JSON string containing a valid, and not already used, UUID. For batch submissions, this value must be a JSON array of valid UUIDs.

If not supplied, the behavior is as today: Cromwell will generate a random workflow ID for every workflow submitted. 

### Bug Fixes

* Fixed a bug on Google Pipelines API backends where missing optional output files (`File?`) were not correctly detected by Cromwell and caused invalid call cache entries to be written.

## 73 Release Notes

### Workflow Restart Performance Improvements

Cromwell now allows for improved performance restarting large workflows through the use of a separate rate limiter for restart checks than the rate limiter used for starting new jobs.
The restart check rate limiter is pre-configured in Cromwell's bundled [reference.conf](https://github.com/broadinstitute/cromwell/blob/develop/core/src/main/resources/reference.conf); see the `job-restart-check-rate-control` stanza in that file for explanations of the various parameters if adjustments are desired.

## 71 Release Notes

### Bug Fixes

* Fixed an issue handling data in Google Cloud Storage buckets with requester pays enabled that could sometimes cause I/O to fail.

## 70 Release Notes

### CWL security fix [#6510](https://github.com/broadinstitute/cromwell/pull/6510)

Fixed an issue that could allow submission of an untrusted CWL file to initiate remote code execution. The vector was improper deserialization of the YAML source file.

CWL execution is enabled by default unless a `CWL` [stanza](https://github.com/broadinstitute/cromwell/blob/develop/core/src/main/resources/reference.conf#L460-L482) is present in the configuration that specifies `enabled: false`. Cromwell instances with CWL disabled were not affected. Consequently, users who wish to mitigate the vulnerability without upgrading Cromwell may do so via this config change.

- Thank you to [Bruno P. Kinoshita](https://github.com/kinow) who first found the issue in a different CWL project ([CVE-2021-41110](https://github.com/common-workflow-language/cwlviewer/security/advisories/GHSA-7g7j-f5g3-fqp7)) and [Michael R. Crusoe](https://github.com/mr-c) who suggested we investigate ours.

## 68 Release Notes

### Virtual Private Cloud

Previous Cromwell versions allowed PAPIV2 jobs to run on a specific subnetwork inside a private network by adding the
information to Google Cloud project labels.

Cromwell now allows PAPIV2 jobs to run on a specific subnetwork inside a private network by adding the network and
subnetwork name directly inside the `virtual-private-cloud` backend configuration. More info
[here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

## 67 Release Notes

### Configuration updates for improved scaling

Some configuration changes were introduced in Cromwell 67 to support improved scaling. See Cromwell's `reference.conf` for details on new parameters.

* I/O throttling moved from `io` to its own `io.throttle` stanza; config updates may be required if these values are currently being overridden in local deployments.

* The default `system.job-rate-control` has been changed from 50 per second to 20 per 10 seconds.

* New configuration parameters have been introduced for values which were previously hardcoded constants:
  * `system.file-hash-batch-size`, value updated from `100` to `50`.
  * `io.gcs.max-batch-size`, value stays the same at `100`.
  * `io.gcs.max-batch-duration`, value stays the same at `5 seconds`.

* New configuration parameters which should not require updating:
  * `io.command-backpressure-staleness`
  * `io.backpressure-extension-log-threshold`
  * `load-control.io-normal-window-minimum`
  * `load-control.io-normal-window-maximum`

* `io.nio.parallelism` was previously misspelled in `reference.conf` but not in Cromwell's configuration reading code. Only correct spellings of this configuration key had or will have effect.

## 66 Release Notes

### Google Artifact Registry Support
Cromwell now supports call caching when using Docker images hosted on
[Google Artifact Registry](https://cloud.google.com/artifact-registry).

### Google Image Repository Hashing Updates
The previously documented `docker.hash-lookup.gcr` configuration has been renamed to `docker.hash-lookup.google` and
now applies to both Google Container Registry (GCR) and Google Artifact Registry (GAR) repositories.
Support for the `docker.hash-lookup.gcr-api-queries-per-100-seconds` configuration key has been formally discontinued
and a bug preventing correct handling of `docker.hash-lookup...throttle` configuration has been fixed.
Please see Cromwell's bundled
[`reference.conf`](https://github.com/broadinstitute/cromwell/blob/develop/core/src/main/resources/reference.conf)
for more details.

## 65 Release Notes

* An additional set of metrics relating to metadata age were added.

### AMD Rome support on PAPI v2
On the PAPI v2 backends "AMD Rome" is now supported as a CPU platform. More details can be found
[here](https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#cpuplatform).

## 64 Release Notes

### Intel Cascade Lake support on PAPI v2

On the PAPI v2 backends "Intel Cascade Lake" is now supported as a CPU platform. More details can be found
[here](https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#cpuplatform).

## 63 Release Notes

### Removed refresh token authentication mode

Google Pipelines API v1 supported authentication with refresh tokens, while v2 of the API does not.

Now that v1 has been discontinued and shut down, this version of Cromwell removes support for refresh tokens.

## 62 Release Notes

### Downloading Access URLs

Added experimental support to download data during Google [Cloud Life Sciences](https://cloud.google.com/life-sciences)
jobs using [DRS
AccessURLs](https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.1.0/docs/#_accessurl).

## 61 Release Notes

### No labels update for Archived workflows

If **- and ONLY if -** you have metadata archiving turned on, then for a workflow whose metadata has been archived by Cromwell 
according to the lifecycle policy, Cromwell will no longer add new labels or update existing labels for this workflow 
coming through PATCH `/labels` endpoint.

## 60 Release Notes

### Java 11

As of this version, a distribution of Java 11 is required to run Cromwell. Cromwell is developed, tested, and
containerized using [AdoptOpenJDK 11 HotSpot](https://adoptopenjdk.net/).

### Hybrid metadata storage ("carboniting") removed

Carboniting functionality has been removed from Cromwell. 
There will be no effect for customers who store metadata permanently in the relational database (most common),
and there will also be no effect for customers who use the in-memory database.

Breaking change only for customers who explicitly enabled `carbonite-metadata-service` in their configuration to split
metadata storage between a relational database and Google Cloud Storage. If you had previously enabled carboniting and 
deletion, any workflows marked as `ArchivedAndPurged` in your database will no longer be accessible via the Cromwell metadata API.

## 59 Release Notes

### Bug Fixes

* Fixed a pair of bugs that could cause workflows to fail unexpectedly with the errors "413 Request Entity Too Large"
  and "java.net.SocketTimeoutException: Read timed out" when accessing Google Cloud Storage.

## 58 Release Notes

Internal CI-related changes only.

## 57 Release Notes

### Breaking configuration change to reference disk support on PAPI v2

Beginning with Cromwell 57, reference disk manifests are now specified completely within Cromwell configuration
rather than through a level of indirection to a manifest file stored in GCS. More details can be found
[here](https://cromwell.readthedocs.io/en/develop/backends/Google#reference-disk-support).

## 56 Release Notes

### Retry with More Memory as workflow option

The experimental memory retry feature gains per-workflow customization and includes breaking changes:
* The per-backend configuration key `<backend>.config.memory-retry.error-keys` has been removed and replaced 
with global key `system.memory-retry-error-keys`
* The per-backend configuration key `<backend>.config.memory-retry.multiplier` has been replaced with **workflow option** 
`memory_retry_multiplier`

More details can be found [here](https://cromwell.readthedocs.io/en/develop/wf_options/Overview.md#retry-with-more-memory-multiplier).

### Bug Fixes

* Fixed a bug that caused Cromwell to mark workflows as failed after a single `500`, `503`, or `504` error from Google Cloud Storage.
  * Cromwell will now retry these errors as designed.
  * The default retry count is `5` and may be customized with `system.io.number-of-attempts`. 

## 55 Release Notes

### Apple Silicon support statement

Users with access to the new Mac hardware should review [important information provided here](https://cromwell.readthedocs.io/en/stable/Releases).

### Bug Fixes

* Fixed a bug that prevented `read_json()` from working with arrays and primitives. The function now works as expected for all valid JSON data inputs. 
More information on JSON Type to WDL Type conversion can be found [here](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#mixed-read_jsonstringfile).

* Now retries HTTP 408 responses as well as HTTP 429 responses during DOS/DRS resolution requests.

* Fixed a bug that prevented the call caching diff endpoint from working with scatters in workflows with archived metadata.

### New Features

#### Reference disk support on PAPI v2

Cromwell now offers support for the use of reference disks on the PAPI v2 backend as an alternative to localizing
reference inputs. More details [here](https://cromwell.readthedocs.io/en/develop/backends/Google#reference-disk-support).

#### Docker image cache support on PAPI v2 lifesciences beta

Cromwell now offers support for the use of Docker image caches on the PAPI v2 lifesciences beta backend. More details [here](https://cromwell.readthedocs.io/en/develop/backends/Google#docker-image-cache-support).

#### Preemptible Recovery via Checkpointing

* Cromwell can now help tasks recover from preemption by allowing them to specify a 'checkpoint' file which will be restored
to the worker VM on the next attempt if the task is interrupted. More details [here](https://cromwell.readthedocs.io/en/develop/optimizations/CheckpointFiles)

## 54 Release Notes

### Bug Fixes

* Fixed a bug that prevented `write_json()` from working with arrays and primitives. The function now works as expected for `Boolean`, `String`, `Integer`, `Float`,
 `Pair[_, _]`, `Object`, `Map[_, _]` and `Array[_]` (including array of objects) type inputs. More information on WDL Type to JSON Type 
 conversion can be found [here](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#mixed-read_jsonstringfile).

### Spark backend support removal

Spark backend was not widely used and it was decided to remove it from the codebase in order to narrow the scope of Cromwell code. 

### Improved DRS Localizer logging

Error logging while localizing a DRS URI should now be more clear especially when there is a Requester Pays bucket involved.

### Per-backend hog factors
Cromwell now allows overriding system-level log factors on back-end level. First, Cromwell will try to use hog-factor 
defined in the backend config, and if it is not defined, it will default to using system-wide hog factor.
```conf
backend {
  providers {
    PAPIv2 {
      config {
        hog-factor: 2
      }
    }
  }
}
```
For more information about hog factors please see [this page](https://cromwell.readthedocs.io/en/develop/cromwell_features/HogFactors/).

### `martha_v2` Support Removed

Cromwell now only supports resolving DOS or DRS URIs through [Martha](https://github.com/broadinstitute/martha)'s most
recent metadata endpoint `martha_v3`, dropping support for Martha's previous metadata endpoint `martha_v2`. To switch to
the new version of Martha's metadata endpoint, update the `martha.url` found in the [filesystems
config](https://cromwell.readthedocs.io/en/stable/filesystems/Filesystems/#overview) to point to `/martha_v3`. More
information on Martha's `martha_v3` request and response schema can be found
[here](https://github.com/broadinstitute/martha#martha-v3).

### DOS/DRS `localization_optional` Support

When running on a backend that supports `localization_optional: true` any DOS or DRS `File` values in the generated
command line will be substituted with the `gsUri` returned from Martha's `martha_v3` endpoint. More information on
`localization_optional` can be found [here](https://cromwell.readthedocs.io/en/stable/optimizations/FileLocalization/).

### DOS/DRS metadata retrieval retried by default

Attempts to retrieve DOS/DRS metadata from Martha will be retried by default. More information can be found
[here](https://cromwell.readthedocs.io/en/stable/optimizations/FileLocalization/).

## 53 Release Notes

### Martha v3 Support

Cromwell now supports resolving DRS URIs through Martha v3 (in addition to Martha v2). To switch to the new version of Martha, update the `martha.url` found in the [filesystems config](https://cromwell.readthedocs.io/en/stable/filesystems/Filesystems/#overview) to
point to `/martha_v3`. More information on Martha v3 request and response schema can be found [here](https://github.com/broadinstitute/martha#martha-v3).

### Support for custom entrypoints on Docker images

Cromwell can now support docker images which have custom entrypoints in the PAPIv2 alpha and beta backends.

### Alpha support for WDL optional outputs on PAPI v2

* Alpha support for WDL optional output files on the PAPI v2 backend has been added, please see the
[documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Google#alpha-support-for-wdl-optional-outputs-on-papi-v2)
for known limitations.

### Monitoring Image Script

* Cromwell now supports an optional `monitoring_image_script` workflow option in addition to the existing
`monitoring_script` and `monitoring_image` options. For more information see the [Google Pipelines API Workflow Options
documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Google#google-pipelines-api-workflow-options).

## 52 Release Notes

### Documentation

Information on how to properly use the Singularity cache with Cromwell is now
provided in the [Cromwell Singularity documentation](
https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#singularity).

### Google library upgrade [(#5565)](https://github.com/broadinstitute/cromwell/pull/5565)

All previous versions of Cromwell shipped with Google Cloud Storage (GCS) libraries that are now deprecated and will [stop working in August 2020](https://developers.googleblog.com/2018/03/discontinuing-support-for-json-rpc-and.html). This release adopts updated libraries to ensure uninterrupted operation. The only user action required is upgrading Cromwell.   

### Bug fixes

* Fixed a bug that required Cromwell to be restarted in order to pick up DNS changes.
    * By default, the JVM caches DNS records with a TTL of infinity.
    * Cromwell now configures its JVM with a 3-minute TTL. This value can be customized by setting `system.dns-cache-ttl`.  
* Clarified an error message that Cromwell emits when the compute backend terminates a job of its own volition (as opposed to termination in response to an abort request from Cromwell)
    * Previously, the error read `The job was aborted from outside Cromwell`
    * The new error reads `The compute backend terminated the job. If this termination is unexpected, examine likely causes such as preemption, running out of disk or memory on the compute instance, or exceeding the backend's maximum job duration.` 

## 51 Release Notes

### Changes and Warnings

The configuration format for call cache blacklisting has been updated, please see the [call caching documentation](
https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching) for details.

### Bug fixes

* Fixed a bug where the `size(...)` function did not work correctly on files 
  from a shared filesystem if `size(...)` was called in the input section on a 
  relative path.
+ Fixed a bug where the `use_relative_output_paths` option would not preserve intermediate folders.

### New functionality

#### Call caching blacklisting improvements

Cromwell previously supported blacklisting GCS buckets containing cache hits which could not be copied for permissions 
reasons. Cromwell now adds support for blacklisting individual cache hits which could not be copied for any reason,
as well as grouping blacklist caches according to a workflow option key. More information available in the [
call caching documentation]( https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching). 

#### new xxh64 and fingerprint strategies for call caching

Existing call cache strategies `path` and `path+modtime` don't work when using docker on shared filesystems 
(SFS backend, i.e. not in cloud storage). The `file` (md5sum) strategy works, but uses a lot of resources.
Two faster strategies have been added for this use case: `xxh64` and 
`fingerprint`. `xxh64` is a lightweight hashing algorithm, `fingerprint` is a strategy designed to be very 
lightweight. Read more about it in the [call caching documentation](
https://cromwell.readthedocs.io/en/stable/Configuring/#call-caching).

## 50 Release Notes

### Changes and Warnings

#### Metadata Archival Config Change

**Note:** Unless you have already opted-in to GCS-archival of metadata during its development, this change will not affect you.
Cromwell's metadata archival configuration has changed in a backwards incompatible way to increase consistency,
please see
[the updated documentation](https://cromwell.readthedocs.io/en/stable/Configuring#hybrid-metadata-storage-classic-carbonite) for details.

## 49 Release Notes

### Changes and Warnings

#### Job store database refactoring

The primary keys of Cromwell's job store tables have been refactored to use a `BIGINT` datatype in place of the previous
`INT` datatype. Cromwell will not be usable during the time the Liquibase migration for this refactor is running.
In the Google Cloud SQL with SSD environment this migration runs at a rate of approximately 40,000 `JOB_STORE_SIMPLETON_ENTRY`
rows per second. In deployments with millions or billions of `JOB_STORE_SIMPLETON_ENTRY` rows the migration may require
a significant amount of downtime so please plan accordingly. The following SQL could be used to estimate the number of
rows in this table:

```
SELECT table_rows FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'cromwell' AND table_name = 'JOB_STORE_SIMPLETON_ENTRY';
```

#### Execution Directory Layout (cache copies)

When an attempt to copy a cache result is made, you'll now see a `cacheCopy` directory in the call root directory. 
This prevents them clashing with the files staged to the same directory for attempt 1 if the cache copy fails (see also: Bug Fixes).

The directory layout used to be:

```
[...]/callRoot/
  - script [from the cache copy attempt, or for execution attempt 1 if the cache copy fails]
  - stdout [from the cache copy attempt, or for execution attempt 1 if the cache copy fails]
  - output.file [from the cache copy attempt, or for execution attempt 1 if the cache copy fails]
  - attempt-2/ [if attempt 1 fails]
    - script
    - stdout
    - output.file
```

but is now:

```
[...]/callRoot/
  - cacheCopy/
    - script
    - stdout
    - output.file
  - script [for attempt 1 if the cache copy fails]
  - stdout [for attempt 1 if the cache copy fails]
  - output.file [for attempt 1 if the cache copy fails]
  - attempt-2/ [if attempt 1 fails]
    - script
    - stdout
    - output.file
```

### New Functionality

#### Disable call-caching for tasks

It is now possible to indicate in a workflow that a task should not be call-cached. See details 
[here](https://cromwell.readthedocs.io/en/stable/optimizations/VolatileTasks).

#### Delete Intermediate Outputs on PapiV2

* **Experimental:** When a new workflow option `delete_intermediate_output_files` is submitted with the workflow,
intermediate `File` objects will be deleted when the workflow completes. See the [Google Pipelines API Workflow Options
documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Google#google-pipelines-api-workflow-options)
for more information.

#### Metadata Archival Support

Cromwell 49 now offers the option to archive metadata to GCS and remove the equivalent metadata from relational
database storage. Please see 
[the documentation](https://cromwell.readthedocs.io/en/stable/Configuring#hybrid-metadata-storage-classic-carbonite) for more details. 

#### Adding support for Google Cloud Life Sciences v2beta
Cromwell now supports running workflows using Google Cloud Life Sciences v2beta API in addition to Google Cloud Genomics v2alpha1. 
More information about migration to the new API from v2alpha1 
[here](https://cromwell.readthedocs.io/en/stable/backends/Google#migration-from-google-cloud-genomics-v2alpha1-to-google-cloud-life-sciences-v2beta). 
* **Note** Google Cloud Life Sciences is the new name for newer versions of Google Cloud Genomics.
* **Note** Support for Google Cloud Genomics v2alpha1 will be removed in a future version of Cromwell. Advance notice will be provided.

### New Docs

#### Installation methods

Links to the conda package and docker container are now available in 
[the install documentation](https://cromwell.readthedocs.io/en/stable/Getting/).


### Bug Fixes

+ Fix a bug where zip files with directories could not be imported. 
  For example a zip with `a.wdl` and `b.wdl` could be imported but one with `sub_workflows/a.wdl` 
  and `imports/b.wdl` could not.
+ Fix a bug which sometimes allowed execution scripts copied by a failed cache-copy to be run instead
  of the attempt-1 script for a live job execution. 
  
## 48 Release Notes

### Womtool Graph for WDL 1.0

The `womtool graph` command now supports WDL 1.0 workflows. 
* **Note:** Generated graphs - including in WDL draft 2 - may look slightly different than they did in version 47.

### Documentation

+ Documented the use of a HSQLDB file-based database so users can try call-caching without needing a database server.
  Please checkout [the database documentation](https://cromwell.readthedocs.io/en/stable/Configuring#database).

## 47 Release Notes

### Retry with more memory on Papiv2 [(#5180)](https://github.com/broadinstitute/cromwell/pull/5180)

Cromwell now allows user defined retries. With `memory-retry` config you can specify an array of strings which when encountered in the `stderr` 
file by Cromwell, allows the task to be retried with multiplier factor mentioned in the config. More information [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### GCS Parallel Composite Upload Support

Cromwell 47 now supports GCS parallel composite uploads which can greatly improve delocalization performance.
This feature is turned off by default, it can be turned on by either a backend-level configuration setting or
on a per-workflow basis with workflow options. More details [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### Papi V2 Localization Using GCR [(#5200)](https://github.com/broadinstitute/cromwell/pull/5200)

The Docker image for the Google Cloud SDK was previously only [published on Docker
Hub](https://hub.docker.com/r/google/cloud-sdk). Now that the image is [publicly hosted in
GCR](http://gcr.io/google.com/cloudsdktool/cloud-sdk), Papi V2 jobs will localize inputs and delocalize outputs using
the GCR image.

## 46 Release Notes

### Nvidia GPU Driver Update

The default driver for Nvidia GPU's on Google Cloud has been updated from `390` to `418.87.00`.  A user may override this option at anytime by providing the `nvidiaDriverVersion` runtime attribute.  See the [Runtime Attribute description for GPUs](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#runtime-attribute-descriptions) for detailed information.

### Enhanced "error code 10" handling in PAPIv2

On Google Pipelines API v2, a worker VM that is preempted may emit a generic error message like
```
PAPI error code 10. The assigned worker has failed to complete the operation
```
instead of a preemption-specific message like
```
PAPI error code 14. Task was preempted for the 2nd time.
```
Cromwell 44 introduced special handling that detects both preemption indicators and re-runs the job consistent with the `preemptible` setting.

Cromwell 46 enhances this handling in response to user reports of possible continued issues.   

## 45 Release Notes

### Improved input and output transfer performance on PAPI v2

Cromwell now requires only a single PAPI "action" each for the entire localization or delocalization process, rather than two per file or directory.
This greatly increases execution speed for jobs with large numbers of input or output files.
In testing, total execution time for a call with 800 inputs improved from more than 70 minutes to less than 20 minutes.

### List dependencies flag in Womtool Command Line [(#5098)](https://github.com/broadinstitute/cromwell/pull/5098)

Womtool now outputs the list of files referenced in import statements using `-l` flag for `validate` command.
More info [here](https://cromwell.readthedocs.io/en/stable/WOMtool/)

### BCS backend new Features support

#### New docker registry
Alibaba Cloud Container Registry is now supported for the `docker` runtime attribute, and the previous `dockerTag` 
runtime attribute continues to be available for Alibaba Cloud OSS Registry.
#### Call caching
Cromwell now supports Call caching when using the BCS backend.
#### Workflow output glob
Globs can be used to define outputs for BCS backend.
#### NAS mount
Alibaba Cloud NAS is now supported for the `mounts` runtime attribute.

### Call Caching Failure Messages [(#5095)](https://github.com/broadinstitute/cromwell/pull/5095)

Call cache failures are no longer sent to the workflow metadata. Instead a limited number of call cache failure messages
will be sent to the workflow log. See [the Cromwell call caching
documentation](https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching/) for more information on call
cache failure logging.

## 44 Release Notes

### Improved PAPI v2 Preemptible VM Support

In some cases PAPI v2 will report the preemption of a VM in a way that differs from PAPI v1. This novel means of reporting
preemption was not recognized by Cromwell's PAPI v2 backend and would result in preemptions being miscategorized as call failures.
Cromwell's PAPI v2 backend will now handle this type of preemption.

## 43 Release Notes

### Virtual Private Cloud with Subnetworks

Cromwell now allows PAPIV2 jobs to run on a specific subnetwork inside a private network by adding the subnetwork key 
`subnetwork-label-key` inside `virtual-private-cloud` in backend configuration. More info [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### Call caching database refactoring

Cromwell's `CALL_CACHING_HASH_ENTRY` primary key has been refactored to use a `BIGINT` datatype in place of the previous
`INT` datatype. Cromwell will not be usable during the time the Liquibase migration for this refactor is running.
In the Google Cloud SQL with SSD environment this migration runs at a rate of approximately 100,000 `CALL_CACHING_HASH_ENTRY`
rows per second. In deployments with millions or billions of `CALL_CACHING_HASH_ENTRY` rows the migration may require  
a significant amount of downtime so please plan accordingly. The following SQL could be used to estimate the number of
rows in this table:

```
select max(CALL_CACHING_HASH_ENTRY_ID) from CALL_CACHING_HASH_ENTRY
```

### Stackdriver Instrumentation

Cromwell now supports sending metrics to [Google's Stackdriver API](https://cloud.google.com/monitoring/api/v3/). 
Learn more on how to configure [here](https://cromwell.readthedocs.io/en/stable/developers/Instrumentation/).

### BigQuery in PAPI

Cromwell now allows a user to specify BigQuery jobs when using the PAPIv2 backend

### Configuration Changes

#### StatsD Instrumentation

There is a small change in StatsD's configuration path. Originally, the path to the config was `services.Instrumentation.config.statsd`
which now has been updated to `services.Instrumentation.config`. More info on its configuration can be found
[here](https://cromwell.readthedocs.io/en/stable/developers/Instrumentation/).

#### cached-copy

A new experimental feature, the `cached-copy` localization strategy is available for the shared filesystem. 
More information can be found in the [documentation on localization](https://cromwell.readthedocs.io/en/stable/backends/HPC).

#### Yaml node limits

Yaml parsing now checks for cycles, and limits the maximum number of parsed nodes to a configurable value. It also
limits the nesting depth of sequences and mappings. See [the documentation on configuring
YAML](https://cromwell.readthedocs.io/en/stable/Configuring/#yaml) for more information.

### API Changes

#### Workflow Metadata

* It is now possible to use `includeKey` and `excludeKey` at the same time. If so, the metadata key must match the `includeKey` **and not** match the `excludeKey` to be included.
* It is now possible to use "`calls`" as one of your `excludeKey`s, to request that only workflow metadata gets returned.

### PostgreSQL support

Cromwell now supports PostgreSQL (version 9.6 or higher, with the Large Object
extension installed) as a database backend.
See [here](https://cromwell.readthedocs.io/en/stable/Configuring/#database) for
instructions for configuring the database connection.

## 42 Release Notes

### Womtool endpoint

The `/describe` endpoint now differentiates between an invalid workflow and a valid workflow with invalid inputs.

Specifically, the new `validWorkflow` key indicates whether the workflow file is valid by itself. If inputs are provided, they are not considered when calculating this field; if inputs are not provided, the value is identical to `valid`.

### Configuration Changes

 *  Virtual private networks can now be configured. See the section below for details.
 
#### Batch Request Timeouts

The timeout on Cromwell's requests to PAPIv2 can now be configured. See the sample PAPIv2.conf for more documentation:

```conf
backend {
  providers {
    PAPIv2 {
      config { 
        batch-requests {
          timeouts {
            read = 10 seconds
            connect = 10 seconds
          }
        }
      }
    }
  }
}
```

### Virtual Private Networks

Cromwell now allows PAPIV2 jobs to run on a private network by adding the network name inside `virtual-private-cloud` in backend configuration.
More info [here](https://cromwell.readthedocs.io/en/stable/backends/Google/).

### AWS Backend

Now includes background job status polling to hopefully reduce the incidence of 'HTTP 429' errors for large workflows.

## 41 Release Notes

### Workflow Options

* It is now possible to supply custom `google-labels` in [workflow options](https://cromwell.readthedocs.io/en/stable/wf_options/Google/).

### AWS backend

It is now possible to use WDL disk attributes with the following formats on AWS.
```
disks: "local-disk 20 SSD"
```
```
disks: "/some/mnt 20 SSD"
```
Because Cromwell's AWS backend auto-sizes disks, the size specification is simply discarded.

### Time Formatting

In previous versions of Cromwell, times were converted to strings using
[the default Java formatter](https://docs.oracle.com/javase/8/docs/api/java/time/OffsetDateTime.html#toString--) which
generates a variety of ISO-8601 formats. String conversions also retained whatever server time zone generated that
specific time instance.

Going forward, times stored in Cromwell metadata, and later returned via the HTTP endpoint, are now converted to UTC
then formatted with exactly three digits of milliseconds.

For example:
- `2017-01-19T12:34:56-04:00` will now be formatted as
- `2017-01-19T16:34:56.000Z`

This change only affects newly formatted dates. Older dates already formatted and stored by previous versions of
Cromwell will not be updated however they will still return a
[valid ISO-8601 format](https://en.wikipedia.org/wiki/ISO_8601). The older format may be in various non-UTC time zones,
and may or may not include microseconds or even nanoseconds, for example `2017-01-19T12:34:56.123456789-04:00`.

### Config Changes

#### Heartbeat failure shutdown

When a Cromwell instance is unable to write heartbeats for some period of time it will automatically shut down. For more
information see the docs on [configuring Workflow Hearbeats](https://cromwell.readthedocs.io/en/stable/Configuring/).

NOTE: In the remote chance that the `system.workflow-heartbeats.ttl` has been configured to be less than `5 minutes`
then the new configuration value `system.workflow-heartbeats.write-failure-shutdown-duration` must also be explicitly
set less than the `ttl`.

#### nVidia Driver Attribute Change

The runtime attribute `nvidia-driver-version` was previously allowed only as a default runtime attribute in configuration.
Because WDL does not allow attribute names to contain `-` characters, this has been changed to `nvidiaDriverVersion`.
This field is now accepted within WDL files as well as within the configuration file.

#### Logging long running jobs

All backends can now emit slow job warnings after a configurable time running. 
NB This example shows how to configure this setting for the PAPIv2 backend:
```conf
# Emit a warning if jobs last longer than this amount of time. This might indicate that something got stuck.
backend {
  providers {
    PAPIv2 {
      config { 
        slow-job-warning-time: 24 hours
      }
    }
  }
}
```

### Runtime Attributes

#### GPU Attributes

* The `gpuType` attribute is no longer validated against a whitelist at workflow submission time. Instead, validation now happens at runtime. This allows any valid accelerator to be used.
* The `nvidiaDriverVersion` attribute is now available in WDL `runtime` sections. The default continues to be `390.46` which applies if and only if GPUs are being used.
* A default `gpuType` ("nvidia-tesla-k80") will now be applied if `gpuCount` is specified but `gpuType` is not.
* Similarly, a default `gpuCount` (1) will be applied if `gpuType` is specified but `cpuCount` is not. 

### Bug fixes

#### Better validation of workflow heartbeats

An error will be thrown on startup when the `system.workflow-heartbeats.heartbeat-interval` is not less than the
`system.workflow-heartbeats.ttl`.


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
### Workflow options changes

A new workflow option is added. If the `final_workflow_outputs_dir` is set 
`use_relative_output_paths` can be used. When set to `true` this will copy 
all the outputs relative to their execution directory. 
my_final_workflow_outputs_dir/~~MyWorkflow/af76876d8-6e8768fa/call-MyTask/execution/~~output_of_interest.
More information can be found in [the workflow options documentation](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/#output-copying).

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

| PAPI Version  |                                 actor-factory                                |
|---------------|:----------------------------------------------------------------------------:|
|      V1       | cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory |
|      V2alpha1 | cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory |
|      V2beta   | cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory   |

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
[cromwell.examples.conf](https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf) for more
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
