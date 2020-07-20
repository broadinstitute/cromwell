**TES Backend**

The TES backend submits jobs to a server that complies with the protocol described by the [GA4GH schema](https://github.com/ga4gh/task-execution-schemas).

This backend creates three files in the `<call_dir>`:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `stdout` - The standard output of the process
* `stderr` - The standard error of the process

The `script` file contains:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

`<container_call_root>` would be equal to the runtime attribute `dockerWorkingDir`  or `/cromwell-executions/<workflow_uuid>/call-<call_name>/execution` if this attribute is not supplied.

**Configuring**

Configuring the TES backend is straightforward; one must only provide the TES API endpoint for the service. 

```hocon
backend {
  default = "TES"
  providers {
    TES {
      actor-factory = "cromwell.backend.impl.tes.TesBackendLifecycleActorFactory"
      config {
        endpoint = "https://<some-url>/v1/tasks"
        root = "cromwell-executions"
        dockerRoot = "/cromwell-executions"
        concurrent-job-limit = 1000
      }
    }
  }
}
```

**Supported File Systems**  

Currently this backend only works with files on a Local or Shared File System. 

**Docker**

This backend supports the following optional [Runtime Attributes](../RuntimeAttributes) and [Workflow Options](../wf_options/Overview/) for working with Docker:

* `docker`: Docker image to use such as "Ubuntu".
* `dockerWorkingDir`: defines the working directory in the container.

**CPU, Memory and Disk** 

This backend supports CPU, memory and disk size configuration through the use of the following [Runtime Attributes](../RuntimeAttributes) and [Workflow Options](../wf_options/Overview/):  

* `cpu` defines the amount of CPU to use. 
    * Type: Integer (ex: 4)
* `memory` defines the amount of memory to use. 
    * Type: String (ex: "4 GB" or "4096 MB")
* `disk` defines the amount of disk to use. 
    * Type: String (ex: "1 GB" or "1024 MB")
* `preemptible` defines whether or not to use preemptible VMs. 
    * Type: Boolean (ex: "true" or "false")

If they are not set, the TES backend may use default values.

**TESK**

[TESK](https://github.com/EMBL-EBI-TSI/TESK) is an implementation of the TES interface that uses Kubernetes and FTP.
When running Cromwell with a TESK backend, you will want to customize the way Cromwell process globs, as kubernetes will not work well with hard links in a lot of cases which is the default behavior in Cromwell.
By adding this to the `config` section of the TES backend in Cromwell, Cromwell will use symlinks instead.  

`glob-link-command = "ls -L GLOB_PATTERN 2> /dev/null | xargs -I ? ln -s ? GLOB_DIRECTORY"`
