Call Caching allows Cromwell to detect when a job has been run in the past so it doesn't have to re-compute results.  Cromwell searches the cache of previously run jobs for a one that has the exact same command and exact same inputs.  If a previously run job is found in the cache, Cromwell will **copy the results** of the previous job instead of re-running it.

Cromwell's call cache is maintained in its database.  For best mileage with call caching, configure Cromwell to [point to a MySQL database](#database) instead of the default in-memory database.  This way any invocation of Cromwell (either with `run` or `server` subcommands) will be able to utilize results from all calls that are in that database.

**Call Caching is disabled by default.**  Once enabled, Cromwell will search the call cache for every `call` statement invocation, assuming `read_from_cache` is enabled (see below):

* If there was no cache hit, the `call` will be executed as normal.  Once finished it will add itself to the cache, assuming `read_from_cache` is enabled (see below)
* If there was a cache hit, outputs are **copied** from the cached job to the new job's output directory

> **Note:** If call caching is enabled, be careful not to change the contents of the output directory for any previously run job.  Doing so might cause cache hits in Cromwell to copy over modified data and Cromwell currently does not check that the contents of the output directory changed.

## Configuring Call Caching
To enable Call Caching, add the following to your Cromwell [configuration](#configuring-cromwell):

```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```

When `call-caching.enabled=true` (default: `false`), Cromwell will be able to to copy results from previously run jobs (when appropriate).
When `invalidate-bad-cache-results=true` (default: `true`), Cromwell will invalidate any cache results which fail to copy during a cache-hit. This is usually desired but might be unwanted if a cache might fail to copy for external reasons, such as a difference in user authentication.

## Call Caching Workflow Options
Cromwell also accepts two [workflow option](#workflow-options) related to call caching:

* If call caching is enabled, but one wishes to run a workflow but not add any of the calls into the call cache when they finish, the `write_to_cache` option can be set to `false`.  This value defaults to `true`.
* If call caching is enabled, but you don't want to check the cache for any `call` invocations, set the option `read_from_cache` to `false`.  This value also defaults to `true`

> **Note:** If call caching is disabled, the workflow options `read_from_cache` and `write_to_cache` will be ignored and the options will be treated as though they were 'false'.

## Docker Tags

Docker tags are a convenient way to point to a version of an image (ubuntu:14.04), or even the latest version (ubuntu:latest).
For that purpose, tags are mutable, meaning that the image they point to can change, while the tag name stays the same.
While this is very convenient in some cases, using mutable, or "floating" tags in tasks affects the reproducibility of a workflow: the same workflow using "ubuntu:latest" run now, and a year, or even a month from now will actually run with different docker images.
This has an even bigger impact when Call Caching is turned on in Cromwell, and could lead to unpredictable behaviors if a tag is updated in the middle of a workflow or even a scatter for example.
Docker provides another way of identifying an image version, using the specific digest of the image. The digest is guaranteed to be different if 2 images have different byte content. For more information see https://docs.docker.com/registry/spec/api/#/content-digests
A docker image with digest can be referenced as follows : **ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950**
The above image refers to a specific image of ubuntu, that does not depend on a floating tag.
A workflow containing this Docker image run now and a year from now will run in the exact same container.

In order to remove unpredictable behaviors, Cromwell takes the following approach regarding floating docker tags.

When Cromwell finds a job ready to be run, it will first look at its docker runtime attribute, and apply the following logic:

* The job doesn't specify a docker image: The job will be dispatched and all call caching settings (read/write) will apply normally.
* The job does specify a docker runtime attribute:
    * The docker image uses a hash: All call caching settings apply normally
    * The docker image uses a floating tag:
        * Cromwell will attempt to look up the hash of the image. Upon success it will pass both the floating tag and this hash value to the backend.
        * All backends currently included with Cromwell will utilize this hash value to run the job.
        * Within a single workflow all floating tags will resolve to the same hash value even if Cromwell is restarted when the workflow is running.
        * If Cromwell fails to lookup the hash (unsupported registry, wrong credentials, ...) it will run the job with the user provided floating tag.
        * The actual Docker image (floating tag or hash) used for the job will be reported in the `dockerImageUsed` attribute of the call metadata.

### Docker Lookup

Cromwell provides 2 methods to lookup a docker hash from a docker tag:

* Local
    In this mode, cromwell will first attempt to find the image on the local machine where it's running using the docker CLI. If the image is present, then its hash will be used.
    If it's not present, cromwell will execute a `docker pull` to try and retrieve it. If this succeeds, the newly retrieved hash will be used. Otherwise the lookup will be considered failed.
    Note that cromwell runs the `docker` CLI the same way a human would. This means two things:
     * The machine Cromwell is running on needs to have docker installed and a docker daemon running.
     * Whichever credentials (and only those) are available on that machine will be available to pull the image.
    
* Remote
    In this mode, cromwell will attempt to retrieve the hash by contacting the remote docker registry where the image is stored. This currently supports Docker Hub and GCR.
    
    Docker registry and access levels supported by Cromwell for docker hash lookup in "remote" mode:
    
    |       |       DockerHub    ||       GCR       ||
    |:-----:|:---------:|:-------:|:------:|:-------:|
    |       |   Public  | Private | Public | Private |
    |  JES  |     X     |    X    |    X   |    X    |
    | Other |     X     |         |    X   |         |

## Local Filesystem Options
When running a job on the Config (Shared Filesystem) backend, Cromwell provides some additional options in the backend's config section:

```
      config {
        ...
        filesystems {
          ...
          local {
            ...
            caching {
              # When copying a cached result, what type of file duplication should occur. Attempted in the order listed below:
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]

              # Possible values: file, path
              # "file" will compute an md5 hash of the file content.
              # "path" will compute an md5 hash of the file path. This strategy will only be effective if the duplication-strategy (above) is set to "soft-link",
              # in order to allow for the original file path to be hashed.
              # Default: file
              hashing-strategy: "file"

              # When true, will check if a sibling file with the same name and the .md5 extension exists, and if it does, use the content of this file as a hash.
              # If false or the md5 does not exist, will proceed with the above-defined hashing strategy.
              # Default: false
              check-sibling-md5: false
            }
          }
        }
      }
```