# Customize tasks

Runtime attributes can be specified in one of two ways:

 1. Within a task you can specify runtime attributes to customize the environment for the call.  
 2. [Default runtime attributes](#default-values) for all tasks can be specified in [Workflow Options](wf_options/Overview.md).

_Task Example_

```
task jes_task {
  command {
    echo "Hello JES!"
  }
  runtime {
    docker: "ubuntu:latest"
    memory: "4G"
    cpu: "3"
    zones: "us-central1-c us-central1-b"
    disks: "/mnt/mnt1 3 SSD, /mnt/mnt2 500 HDD"
  }
}
workflow jes_workflow {
  call jes_task
}
```


## Recognized Runtime attributes and Backends

Cromwell recognizes certain runtime attributes and has the ability to format these for some [Backends](/backends/Backends). See the table below for common attributes that apply to _most_ backends.

| Runtime Attribute                               | Local | Google Cloud | TES       | AWS Batch |            HPC            |
|-------------------------------------------------|:-----:|:------------:|-----------|:---------:|:-------------------------:|
| [`cpu`](#cpu)                                   |       |      ✅       |           |     ✅     |           `cpu`           |
| [`memory`](#memory)                             |       |      ✅       |           |     ✅     | `memory_mb` / `memory_gb` |
| [`disks`](#disks)                               |       |      ✅       | ⚠️ Note 1 | ⚠️ Note 2 |         ℹ️ Note 3         |
| [`disk`](#disk)                                 |       |              | ✅         |           |                           |
| [`docker`](#docker)                             |   ✅   |      ✅       |           |     ✅     |    `docker` ℹ️ Note 3     |
| [`maxRetries`](#maxretries)                     |   ✅   |      ✅       |           |     ✅     |         ℹ️ Note 3         |
| [`continueOnReturnCode`](#continueonreturncode) |   ✅   |      ✅       |           |     ✅     |         ℹ️ Note 3         |
| [`failOnStderr`](#failonstderr)                 |   ✅   |      ✅       |           |     ✅     |         ℹ️ Note 3         |


> **Note 1**
> 
> Partial support. See [TES documentation](/backends/TES) for details. 
 
> **Note 2**
>
> Partial support. See [`disks`](#disks) for details.

> **Note 3**
> 
> The HPC [Shared Filesystem backend](/backends/HPC#shared-filesystem) (SFS) is fully configurable and any number of attributes can be exposed. Cromwell recognizes some of these attributes (`cpu`, `memory` and `docker`) and parses them into the attribute listed in the table which can be used within the HPC backend configuration.


### Google Cloud Specific Attributes
There are a number of additional runtime attributes that apply to the Google Cloud Platform:

- [zones](#zones)
- [preemptible](#preemptible)
- [bootDiskSizeGb](#bootdisksizegb)
- [noAddress](#noaddress)
- [gpuCount, gpuType, and nvidiaDriverVersion](#gpucount-gputype-and-nvidiadriverversion)
- [cpuPlatform](#cpuplatform)
- [useDockerImageCache](#usedockerimagecache)



## Expression support

Runtime attribute values are interpreted as expressions.  This means that it has the ability to express the value of a runtime attribute as a function of one of the task's inputs.  
_For example:_

```
task runtime_test {
  String ubuntu_tag
  Int memory_gb

  command {
    ./my_binary
  }

  runtime {
    docker: "ubuntu:" + ubuntu_tag
    memory: memory_gb + "GB"
  }
}
```

HPC backends may define other configurable runtime attributes beyond the five listed, to find out more visit the [SunGridEngine](/backends/SGE) tutorial.

## Default Values

Default values for runtime attributes can be specified via [Workflow Options](wf_options/overview).  
For example, consider this WDL file:

```wdl
task first {
  command { ... }
}

task second {
  command {...}
  runtime {
    docker: "my_docker_image"
  }
}

workflow w {
  call first
  call second
}
```

And this set of workflow options:

```json
{
  "default_runtime_attributes": {
    "docker": "ubuntu:latest",
    "zones": "us-central1-c us-central1-b"
  }
}
```

Then, these values for `docker` and `zones` will be used for any task that does not explicitly override them in the WDL file. In return, the effective runtime for `task first` is:

```
{
    "docker": "ubuntu:latest",
    "zones": "us-central1-c us-central1-b"
  }
```

And the effective runtime for `task second` is:

```
{
    "docker": "my_docker_image",
    "zones": "us-central1-c us-central1-b"
  }
```

Note how for `task second` the WDL value for `docker` is used instead of the default provided in the workflow options.


## Runtime Attribute Descriptions

### `cpu`

*Default: _1_*

The `cpu` runtime attribute represents the number of cores that a job requires, however each backend may interpret this differently:

- In Google Cloud: this is interpreted as "the minimum number of cores to use."
- In HPCs (SFS): this is configurable, but usually a reservation and/or limit of number of cores.

Example
```
runtime {
  cpu: 2
}
```

### `memory`
*Default: "2G"*

Memory is the amount of RAM that should be allocated to a task, however each backend may interpret this differently:

- Google Cloud: The minimum amount of RAM to use.
- SFS: Configurable, but usually a reservation and/or limit of memory.

The memory size is specified as an amount and units of memory, for example "4G":

```
runtime {
  memory: "4G"
}
```

Within the SFS backend, you can additionally specify `memory_mb` or `memory_gb` as runtime attributes within the configuration. More information can be found [here](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/#specifying-the-runtime-attributes-for-your-hpc-tasks).


### `disks`

This attribute specifies volumes that will be mounted to the VM for your job. These volumes are where you can read and write files that will be used by the commands within your workflow. 


They are specified as a comma separated list of disks. Each disk is further separated as a space separated triplet (e.g. `local-disk 10 SSD`) consisting of:

1. Mount point (absolute path), or `local-disk` to reference the mount point where Google Cloud will localize files and the task's current working directory will be
2. Disk size in GB (rounded to the next 375 GB for LOCAL)
3. Disk type.  One of: "LOCAL", "SSD", or "HDD" ([documentation](https://cloud.google.com/compute/docs/disks/#overview))

All tasks launched on Google Cloud *must* have a `local-disk`.  If one is not specified in the runtime section of the task, then a default of `local-disk 10 SSD` will be used.  The `local-disk` will be mounted to `/cromwell_root`.

For the AWS Batch backend, the disk volume is managed by AWS EBS with autoscaling capabilities.  As such, the Disk size and disk type will be ignored. If provided, the mount point will be verified at runtime.


The Disk type must be one of "LOCAL", "SSD", or "HDD". When set to "LOCAL", the size of the drive is constrained to 375 GB intervals so intermediate values will be rounded up to the next 375 GB. All disks are set to auto-delete after the job completes.

*Example 1: Changing the Localization Disk*

```
runtime {
  disks: "local-disk 100 SSD"
}
```

*Example 2: Mounting an Additional Two Disks*

```
runtime {
  disks: "/mnt/my_mnt 3 SSD, /mnt/my_mnt2 500 HDD"
}
```

### `disk`

Specific to the TES backend, sets the `disk_gb` resource.

```
runtime {
  disk: "25 GB"
}
```

### `docker`

When specified, Cromwell will run your task within the specified Docker image. 

```
runtime {
  docker: "ubuntu:latest"
}
```

- Local: Cromwell will automatically run the docker container.
- SFS: When a docker container exists within a task, the `submit-docker` method is called. See the [Getting started with containers](/tutorials/Containers/) guide for more information.
- GCP: This attribute is mandatory when submitting tasks to Google Cloud.
- AWS Batch: This attribute is mandatory when submitting tasks to AWS Batch.


### `maxRetries`

*Default: _0_*

This retry option is introduced to provide a method for tackling transient job failures. For example, if a task fails due to a timeout from accessing an external service, then this option helps re-run the failed the task without having to re-run the entire workflow. It takes an Int as a value that indicates the maximum number of times Cromwell should retry a failed task. This retry is applied towards jobs that fail while executing the task command. This method only applies to transient job failures and is a feeble attempt to retry a job, that is it cannot be used to increase memory in out-of-memory situations.

If using the Google backend, it's important to note that The `maxRetries` count is independent from the [preemptible](#preemptible) count. For example, the task below can be retried up to 6 times if it's preempted 3 times AND the command execution fails 3 times.

```
runtime {
  preemptible: 3
  maxRetries: 3
}
```

### `continueOnReturnCode`
*Default: _0_*

When each task finishes it returns a code. Normally, a non-zero return code indicates a failure. However you can override this behavior by specifying the `continueOnReturnCode` attribute.

When set to false, any non-zero return code will be considered a failure. When set to true, all return codes will be considered successful.

```
runtime {
  continueOnReturnCode: true
}
```

When set to an integer, or an array of integers, only those integers will be considered as successful return codes.

```
runtime {
  continueOnReturnCode: 1
}
```

```
runtime {
  continueOnReturnCode: [0, 1]
}
```

### `failOnStderr`

*Default: _false_*

Some programs write to the standard error stream when there is an error, but still return a zero exit code. Set `failOnStderr` to true for these tasks, and it will be considered a failure if anything is written to the standard error stream.

```
runtime {
  failOnStderr: true
}
```



### `zones`

The ordered list of zone preference (see [Region and Zones](https://cloud.google.com/compute/docs/zones) documentation for specifics).

*The zones are specified as a space separated list, with no commas:*

```
runtime {
  zones: "us-central1-c us-central1-b"
}
```

Defaults to the configuration setting `genomics.default-zones` in the Google Cloud configuration block, which in turn defaults to using `us-central1-b`.

### `preemptible`

*Default: _0_*

Passed to Google Cloud: "If applicable, preemptible machines may be used for the run."

Take an Int as a value that indicates the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one.  
*eg. With a value of 1, Cromwell will request a preemptible VM, if the VM is preempted, the task will be retried with a non-preemptible VM.*

```
runtime {
  preemptible: 1
}
```

In GCP Batch, preempted jobs can be identified in job metadata (`gcloud batch jobs describe`) by a `statusEvent` with a description that looks like:
```
Job state is set from RUNNING to FAILED for job projects/abc/locations/us-central1/jobs/job-abc.Job
failed due to task failure. Specifically, task with index 0 failed due to the
following task event: "Task state is updated from RUNNING to FAILED on zones/us-central1-b/instances/8675309
due to Spot VM preemption with exit code 50001."
```


### `bootDiskSizeGb`

In addition to working disks, Google Cloud allows specification of a boot disk size. This is the disk where the docker image itself is booted (**not the working directory of your task on the VM**).
Its primary purpose is to ensure that larger docker images can fit on the boot disk.
```
runtime {
  # Yikes, we have a big OS in this docker image! Allow 50GB to hold it:
  bootDiskSizeGb: 50
}
```

Since no `local-disk` entry is specified, Cromwell will automatically add `local-disk 10 SSD` to this list.


### `noAddress`

This runtime attribute adds support to disable assigning external IP addresses to VMs provisioned by the Google backend. If set to true, the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs.

Note well!  You must enable "Private Google Access" for this feature to work. See "How To Setup" below.

For example, the task below will succeed:
```
command {
  echo "hello!"
  
}

runtime {
  docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
  noAddress: true
}
```

The task below will fail for two reasons:
 1. The command is accessing an external service, in this case GitHub.
 2. The docker image is available in DockerHub and not the Google Container Registry. 
```
command {
  git clone https://github.com/broadinstitute/cromwell.git
  
}

runtime {
  docker: "docker.io/alpine/git:latest"
  noAddress: true
}
```

#### How to Setup

Configure your Google network to use "Private Google Access". This will allow your VMs to access Google Services including Google Container Registry, as well as Dockerhub images.

1. Using `gcloud compute networks subnets list`, identify the subnet and region you will be using with Cromwell. If multiple, run the next step for each region and subnet you wish to use.
1. `gcloud compute networks subnets update [SUBNET-NAME] --region [REGION]  --enable-private-ip-google-access`

That's it!  You can now run with `noAddress` runtime attribute and it will work as expected.

### `gpuCount`, `gpuType`, and `nvidiaDriverVersion`

Attach GPUs to the instance when running on the Pipelines API([GPU documentation](https://cloud.google.com/compute/docs/gpus/)).
Make sure to choose a zone for which the type of GPU you want to attach is available.

The types of compute GPU supported are:

* `nvidia-tesla-v100`
* `nvidia-tesla-p100`
* `nvidia-tesla-p4`
* `nvidia-tesla-t4`

On Life Sciences API, the default driver is `418.87.00`. You may specify your own via the `nvidiaDriverVersion` key.  Make sure that driver exists in the `nvidia-drivers-us-public` beforehand, per the [Google Pipelines API documentation](https://cloud.google.com/genomics/reference/rest/Shared.Types/Metadata#VirtualMachine). 

On GCP Batch, `nvidiaDriverVersion` is currently ignored; Batch selects the correct driver version automatically.

```
runtime {
    gpuType: "nvidia-tesla-t4"
    gpuCount: 2
    nvidiaDriverVersion: "418.87.00"
    zones: ["us-central1-c"]
}
```

### `cpuPlatform`

This option is specific to the Google Cloud backend, specifically [this](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform) feature when a certain minimum CPU platform is desired.

A usage example:

```
runtime {
    cpu: 2
    cpuPlatform: "Intel Cascade Lake"
}
```
Note that when this options is specified, make sure the requested CPU platform is [available](https://cloud.google.com/compute/docs/regions-zones/#available) in the `zones` you selected.

The following CPU platforms are currently supported by the Google Cloud backend:
- `Intel Ice Lake`
- `Intel Cascade Lake`
- `Intel Skylake`     
- `Intel Broadwell`   
- `Intel Haswell`     
- `Intel Ivy Bridge`  
- `Intel Sandy Bridge`
- `AMD Rome`

### 'useDockerImageCache'

This option is specific to the Google Cloud backend, moreover it is only supported by Google Life Sciences API starting from version v2 beta.
In order to use this feature Cromwell has to have PAPI v2 backend configured with this feature enabled.  
More information about this feature and it's configuration can be found [in the Google backend section of documentation](backends/Google.md).
