#Backends

A backend is a way to run the commands of your workflow. Cromwell allows for backends conforming to
the Cromwell backend specification to be plugged into the Cromwell engine. Additionally, backends are included with the
Cromwell distribution:

* **[Local](Local)**
* **[HPC](HPC)**, including **[Sun Grid Engine](SGE), [LSF](LSF), [HTCondor](HTcondor) & [SLURM](SLURM)** 
    * Run jobs as subprocesses or via a dispatcher.
    * Supports launching in Docker containers.
    * Use `bash`, `qsub`, and `bsub` to run scripts.
* **[Google Cloud](Google)** 
    * Launch jobs on Google Compute Engine through the Google Genomics Pipelines API.
* **[GA4GH TES](TES)** 
    * Launch jobs on servers that support the GA4GH Task Execution Schema (TES).
* **[Spark](Spark)** 
    * Supports execution of Spark jobs.
* **[Alibaba Cloud](BCS)** 
    * Launch jobs on Alibaba Cloud BatchCompute service.
* **[AWS Batch (beta)](AWS.md)**
    * Use Job Queues on AWS Batch

HPC backends are put under the same umbrella because they all use the same generic configuration that can be specialized to fit the need of a particular technology.

Backends are specified in the `backend.providers` configuration. Each backend has a configuration that looks like:

```hocon
BackendName {
  actor-factory = "FQN of BackendLifecycleActorFactory class"
  config {
    ...
  }
}
```

The structure within the `config` block will vary from one backend to another; it is the backend implementation's responsibility
to be able to interpret its configuration.

The providers section can contain multiple backends which will all be available to Cromwell.

## Backend Job Limits

All backends support limiting the number of concurrent jobs by specifying the following option in the backend's configuration
stanza:

```
backend {
  ...
  providers {
    BackendName {
      actor-factory = ...
      config {
        concurrent-job-limit = 5
```

## Backend Filesystems

Each backend will utilize a filesystem to store the directory structure and results of an executed workflow.
The backend/filesystem pairings are as follows:

* Local, HPC and Spark backend use the [Shared Local Filesystem](HPC/#filesystems).
* Google backend uses the [Google Cloud Storage Filesystem](Google/#google-cloud-storage-filesystem).
* Alibaba Cloud backend uses the OSS Storage FileSystem.

Additional filesystems capabilities can be added depending on the backend.
For instance, an HPC backend can be configured to work with files on Google Cloud Storage. See the [HPC documentation](HPC) for more details.
