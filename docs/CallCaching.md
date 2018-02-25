Call Caching allows Cromwell to detect when a job has been run in the past so that it doesn't have to re-compute results, saving both time and money.  Cromwell searches the cache of previously run jobs for one that has the exact same command and exact same inputs.  If a previously run job is found in the cache, Cromwell will use the results of the previous job instead of re-running it.

Cromwell's call cache is maintained in its database.  In order for call caching to be used on any previously run jobs, it is best to configure Cromwell to [point to a MySQL database](Configuring#database) instead of the default in-memory database.  This way any invocation of Cromwell (either with `run` or `server` subcommands) will be able to utilize results from all calls that are in that database.

**Configuring Call Caching**

*Call Caching is disabled by default.*  Call Caching can be enabled in your Cromwell [Configuration](Configuring#call-caching) and the behavior can be modified via [Workflow Options](wf_options/Overview). If you are adding Workflow options, do not set [`read_from_cache` or `write_to_cache`](wf_options/Overview#call-caching-options) = false, as it will impact the following process.

Once enabled, Cromwell will search the call cache for every `call` statement invocation.

* If there was no cache hit, the `call` will be executed as normal.  Once finished it will add itself to the cache.
* If there was a cache hit, outputs are either **copied from the original cached job to the new job's output directory** or **referenced from the original cached job** depending on the Cromwell [Configuration](Configuring#call-caching) settings.

> **Note:** If call caching is enabled, be careful not to change the contents of the output directory for any previously run job.  Doing so might cause cache hits in Cromwell to copy over modified data and Cromwell currently does not check that the contents of the output directory changed.  Additionally, if any files from a previous job directory are removed, call caching will fail due to missing files.

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
    
    |       |       DockerHub    ||       GCR       ||
    |:-----:|:---------:|:-------:|:------:|:-------:|
    |       |   Public  | Private | Public | Private |
    | Pipelines API  |     X     |    X    |    X   |    X    |
    | Other |     X     |         |    X   |         |
