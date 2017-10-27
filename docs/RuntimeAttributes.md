# Runtime Attributes: Used to customize tasks

They can be specified in one of two ways:

 1. Within a task you can specify runtime attributes to customize the environment for the call.  
 2. [Default runtime attributes](#default-values) for all tasks can be specified in [Workflow Options](wf_options/Overview.md).

 >* Certain [Backends](backends/Backends) only support certain runtime attributes.  See [Backend Support](#backend-support) for a table.

# Task Example

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

# Expression support

Runtime attribute values are interpreted as expressions.  This means that it has the ability to express the value of a runtime attribute as a function of one of the task's inputs.  
For example:

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

SGE and similar backends may define other configurable runtime attributes beyond the five listed. To find more information about SGE, view [Sun GridEngine](backends/SGE).

# Default Values

Default values for runtime attributes can be specified via [Workflow Options](wf_options/WorkflowOptions).  
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
Note how for `task second`:  
*The WDL value for `docker` is used instead of the default provided in the workflow options.*

## `continueOnReturnCode`

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

Defaults to "0".

## `cpu`

Passed to Google Cloud: "The minimum number of cores to use."

Passed to SGE, etc.: Configurable, but usually a reservation and/or limit of number of cores.

```
runtime {
  cpu: 2
}
```

Defaults to "1".

## `disks`

This is currently used by the [Google Cloud backend](backends/Google). You can use this attribute to specify volumes that will be mounted to the VM for your job.  These volumes are where you can read and write files that will be used by the commands within your workflow.

They are specified as a comma separated list of disks. Each disk is further separated as a space separated triplet (e.g. `local-disk 10 SSD`) consisting of:

1. Mount point (absolute path), or `local-disk` to reference the mount point where Google Cloud will localize files and the task's current working directory will be
2. Disk size in GB (ignored for disk type LOCAL)
3. Disk type.  One of: "LOCAL", "SSD", or "HDD" ([documentation](https://cloud.google.com/compute/docs/disks/#overview))

All tasks launched on Google Cloud *must* have a `local-disk`.  If one is not specified in the runtime section of the task, then a default of `local-disk 10 SSD` will be used.  The `local-disk` will be mounted to `/cromwell_root`.

The Disk type must be one of "LOCAL", "SSD", or "HDD". When set to "LOCAL", the size of the drive is automatically provisioned by Google so any size specified in WDL will be ignored. All disks are set to auto-delete after the job completes.

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

## `bootDiskSizeGb`

In addition to working disks, Google Cloud allows specification of a boot disk size. This is the disk where the docker image itself is booted (**not the working directory of your task on the VM**).
Its primary purpose is to ensure that larger docker images can fit on the boot disk.
```
runtime {
  # Yikes, we have a big OS in this docker image! Allow 50GB to hold it:
  bootDiskSizeGb: 50
}
```

Since no `local-disk` entry is specified, Cromwell will automatically add `local-disk 10 SSD` to this list.

## `zones`

The ordered list of zone preference (see [Region and Zones](https://cloud.google.com/compute/docs/zones) documentation for specifics)

*The zones are specified as a space separated list, with no commas:*

```
runtime {
  zones: "us-central1-c us-central1-b"
}
```

Defaults to the configuration setting `genomics.default-zones` in the Google Cloud configuration block which in turn defaults to using `us-central1-b`

## `docker`

When specified, cromwell will run your task within the specified Docker image.

```
runtime {
  docker: "ubuntu:latest"
}
```

*This attribute is mandatory when submitting tasks to Google Cloud. When running on other backends, they default to not running the process within Docker.*

## `failOnStderr`

Some programs write to the standard error stream when there is an error, but still return a zero exit code. Set `failOnStderr` to true for these tasks, and it will be considered a failure if anything is written to the standard error stream.

```
runtime {
  failOnStderr: true
}
```

*Defaults to "false".*  

## `memory`

Passed to Google Cloud: "The minimum amount of RAM to use."

Passed to SGE, etc.: Configurable, but usually a reservation and/or limit of memory.

The memory size is specified as an amount and units of memory, for example "4 G".

```
runtime {
  memory: "4G"
}
```

*Defaults to "2G".*

## `preemptible`

Passed to Google Cloud: "If applicable, preemptible machines may be used for the run."

Take an Int as a value that indicates the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one.  
*eg. With a value of 1, Cromwell will request a preemptible VM, if the VM is preempted, the task will be retried with a non-preemptible VM.*

```
runtime {
  preemptible: 1
}
```

*Defaults to 0.*

# Backend Support

[Backends](backends/Backends) only support certain attributes. See table below:

| Runtime Attribute    | LOCAL |  Google Cloud  |  SGE  |
| -------------------- |:-----:|:-----:|:-----:|
| [continueOnReturnCode](#continueonreturncode) |   x   |   x   |   x   |
| [cpu](#cpu)                  |       |   x   |   x   |
| [disks](#disks)                                |       |   x   |       |
| [zones](#zones)                                |       |   x   |       |
| [docker](#docker)                              |   x   |   x   |   x   |
| [failOnStderr](#failOnStderr)                  |   x   |   x   |   x   |
| [memory](#memory)                              |       |   x   |   x   |
| [preemptible](#preemptible)                    |       |   x   |       |
| [bootDiskSizeGb](#bootdisksizegb)              |       |   x   |       |

- [Shared Filesystem backend](backends/SharedFilesystem) is fully configurable and thus these attributes do not apply universally

TODO:

- Document TES, HTCondor, Spark backend attribute support
