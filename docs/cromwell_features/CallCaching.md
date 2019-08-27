Call Caching allows Cromwell to detect when a job has been run in the past so that it doesn't have to re-compute results, saving both time and money.  Cromwell searches the cache of previously run jobs for one that has the exact same command and exact same inputs.  If a previously run job is found in the cache, Cromwell will use the results of the previous job instead of re-running it.

Cromwell's call cache is maintained in its database.  In order for call caching to be used on any previously run jobs,
it is best to configure Cromwell to [point to a MySQL database](../../Configuring/#database) instead of the default
in-memory database.  This way any invocation of Cromwell (either with `run` or `server` subcommands) will be able to
utilize results from all calls that are in that database.

**Configuring Call Caching**

*Call Caching is disabled by default.*  Call Caching can be enabled in your Cromwell
[Configuration](../../Configuring/#call-caching) and the behavior can be modified via
[Workflow Options](../../wf_options/Overview/). If you are adding Workflow options, do not set
[`read_from_cache` or `write_to_cache`](../../wf_options/Overview/#call-caching-options) = false, as it will impact the
following process.

Once enabled, Cromwell by default will search the call cache for every `call` statement invocation.

* If there was no cache hit, the `call` will be executed as normal.  Once finished it will add itself to the cache.
* If there was a cache hit, outputs are either **copied from the original cached job to the new job's output directory**
or **referenced from the original cached job** depending on the Cromwell
[Configuration](../../Configuring/#call-caching) settings.

> **Note:** If call caching is enabled, be careful not to change the contents of the output directory for any previously run job.  Doing so might cause cache hits in Cromwell to copy over modified data and Cromwell currently does not check that the contents of the output directory changed.  Additionally, if any files from a previous job directory are removed, call caching will fail due to missing files.

***File hash caching***

Cromwell offers the option to cache file hashes within the scope of a root workflow to prevent repeatedly requesting the hashes of the
same files multiple times. File hash caching is off by default and can be turned on with the configuration option `system.file-hash-cache=true`.

***Call cache copy authorization failure prefix blacklisting***

Cromwell has the option to filter call cache hits based on authorization failures copying previous 
call cache hits. In a multi-user environment user A might cache hit to one of user B's results
but that doesn't necessarily mean user A is authorized to read user B's outputs from the filesystem. Call cache blacklisting
allows Cromwell to record on a per-root-workflow level which file path prefixes were involved in cache result copy authorization failures.
If Cromwell sees that the file paths for a candidate cache hit have a blacklisted prefix, Cromwell will quickly 
fail the copy attempt without doing any potentially expensive I/O.

Call cache blacklisting configuration looks like:

```
call-caching {

  enabled = true

  # In a multi-user environment this should be false so unauthorized users don't invalidate results for authorized users. 
  invalidate-bad-cache-results = false

  blacklist-cache {
    # The call caching blacklist cache is off by default. This is used to blacklist cache hit paths based on the
    # prefixes of cache hit paths that Cromwell previously failed to copy for authorization reasons.
    enabled: true
    # Guava cache concurrency.
    concurrency: 10000
    # How long entries in the cache should live from the time of their last access.
    ttl: 20 minutes
    # Maximum number of entries in the cache.
    size: 1000
  }
}
```

Call cache blacklisting could be supported by any backend type though is currently implemented only for the Google Pipelines API (PAPI) backends.
For PAPI backends the bucket is considered the prefix for blacklisting purposes.

***Call cache hit path prefixes***
 
In a multi-user environment where access to job outputs may be restricted among different users, it can be useful to limit
cache hits to those that are more likely to actually be readable for cache hit copies.
Cromwell now supports a `call_cache_hit_path_prefixes` workflow option for this purpose. This is particularly useful in the PAPI backend where the workflow
root can be specified in workflow options via `jes_gcs_root`. The value of `call_cache_hit_path_prefixes` should be an array of strings representing  
prefixes that call cache hit output files should have in order to be considered as a cache hit. Using PAPI as an example and assuming Alice and Bob have
made their data accessible to each other, Alice could submit a workflow with these options:

```
{
  "call_cache_hit_path_prefixes": [ "gs://alice_bucket", "gs://bob_bucket" ]
}
```

With these workflow options Cromwell would only look for cache hits for Alice's jobs in Alice's or Bob's buckets.

As a further optimization the PAPI backend has the concept of "this" bucket on a per-workflow basis, where "this" bucket is
the bucket that contains the current workflow root.
If `call_cache_hit_path_prefixes` is specified in 
workflow options on the PAPI backend, Cromwell will automatically prepend "this" bucket to the call cache hit path prefixes to search.
For example, if Charles specified the same workflow options as in the example above and his workflow root was under `gs://charles_bucket`,
Cromwell would search cache hits in all of the `gs://alice_bucket`, `gs://bob_bucket` and `gs://charles_bucket` buckets without having to specify
`gs://charles_bucket` bucket explicitly in `call_cache_hit_path_prefixes`.

If no `call_cache_hit_path_prefixes` are specified then all matching cache hits will be considered.

***Call cache failure logging***

When Cromwell fails to cache a job from a previous result the reason will be logged. To reduce the verbosity of the logs
only the first three failure reasons will be logged per shard of each job. Cromwell will continue to try copying
previous results for the call, and when no candidates are left Cromwell will run the job on the backend.

**Docker Tags**

Certain Docker tags can impact call caching performance. 
Docker tags are a convenient way to point to a version of an image (`ubuntu:14.04`), or even the latest version (`ubuntu:latest`).
For that purpose, tags are mutable, meaning that the image they point to can change, while the tag name stays the same.
While this is very convenient in some cases, using mutable, or "floating" tags in tasks affects the reproducibility of a workflow. 
If you were to run the same workflow using `ubuntu:latest` now, and again in a year (or even in a month) may run with different docker images.
This has an even bigger impact when Call Caching is turned on in Cromwell, and could lead to unpredictable behaviors if a tag is updated in the middle of a workflow or even a scatter for example.

In order to ensure perfect reproducibility, Docker provides another way of identifying an image version by using the specific digest of the image, which is an immutable identifier. The digest is guaranteed to be different if 2 images have different byte content. For more information see [Docker's api specs](https://docs.docker.com/registry/spec/api/#/content-digests).
A docker image can be referenced using the digest (e.g. `ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950`).
This image refers to a specific image of ubuntu that does not depend on a floating tag.
A workflow containing this Docker image run now and a year from now will run in the exact same container.

However, in order to remove unpredictable behaviors, when Cromwell finds a job ready to be run, it will first look at its docker runtime attribute, and apply the following logic:

* If the job doesn't specify a docker image it will be dispatched and all call caching settings (read/write) will apply normally.
* If the job does specify a docker runtime attribute, then:
    * if the docker image uses a hash, all call caching settings apply normally
    * if the docker image uses a floating tag:
        * Cromwell will attempt to look up the immutable digest of the image with this floating tag. Upon success it will pass both the floating tag and this digest value to the backend.
        * All backends currently included with Cromwell will utilize this digest value to run the job.
        * Within a single workflow, all floating tags within a given workflow will resolve to the same digest value even if Cromwell is restarted when the workflow is running.
        * If Cromwell fails to lookup the digest (for instance an unsupported docker registry, wrong credentials, ...) it will run the job with the user provided floating tag.
        * The actual docker image (floating tag or digest) used for the job will be reported in the `dockerImageUsed` attribute of the call metadata.

**Docker Lookup**

Cromwell provides two methods to lookup a Docker hash from a Docker tag:

* _Local_  
    In this mode, Cromwell will first attempt to find the image on the local machine where it's running using the `docker` CLI. If the image is locally present, then its digest will be used.
    If the image is not present locally, Cromwell will execute a `docker pull` to try and retrieve it. If this succeeds, the newly retrieved digest will be used. Otherwise the lookup will be considered failed.
    Note that Cromwell runs the `docker` CLI the same way a human would. This means two things:
     * The machine Cromwell is running on needs to have Docker installed and a Docker daemon running.
     * The current `docker` CLI credentials on that machine will be used to pull the image.
    
* _Remote_  
    In this mode, Cromwell will attempt to retrieve the hash by contacting the remote docker registry where the image is stored. This currently supports Docker Hub and GCR.
    
    Docker registry and access levels supported by Cromwell for docker digest lookup in "remote" mode:
    
    <!-- Pasted into then regenerated at https://www.tablesgenerator.com/markdown_tables -->

    |               | DockerHub | DockerHub |   GCR  |   GCR   |   ECR  |   ECR   |   ACR  |   ACR   |
    |:-------------:|:---------:|:---------:|:------:|:-------:|:------:|:-------:|:------:|:-------:|
    |               |   Public  |  Private  | Public | Private | Public | Private | Public | Private |
    | Pipelines API |     X     |     X     |    X   |    X    |        |         |        |         |
    |   AWS Batch   |     X     |           |    X   |         |        |         |        |         |
    |      BCS      |           |           |        |         |        |         |        |    X    |
    |     Other     |     X     |           |    X   |         |        |         |        |         |

    <!-- Pasted then regenerated at https://www.tablesgenerator.com/markdown_tables -->

**Runtime Attributes**

As well as call inputs and the command to run, call caching considers the following [runtime
attributes](../../RuntimeAttributes/) of a given task when determining whether to call cache:

* [`ContinueOnReturnCode`](../../RuntimeAttributes/#continueonreturncode)
* [`Docker`](../../RuntimeAttributes/#docker)
* [`FailOnStderr`](../../RuntimeAttributes/#failonstderr)

If any of these attributes have changed from a previous instance of the same task, that instance will not be call-cached
from. Other runtime attributes, including [`memory`](../../RuntimeAttributes/#memory),
[`cpu`](../../RuntimeAttributes/#cpu), and [`disks`](../../RuntimeAttributes/#disks), are not considered by call caching
and therefore may be changed without preventing a cached result from being used.
