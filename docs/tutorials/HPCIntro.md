## Getting started on HPC clusters

### Prerequisites

This tutorial page relies on completing the previous tutorials:

* [Configuration Files](ConfigurationFiles.md)

### Goals

At the end of this tutorial you'll have set up Cromwell to run against your HPC cluster. We'll use SGE as an example but this applies equally to LSF and others.

### Let's get started!

####

#### Telling Cromwell the type of backend

Start by defining your new backend configuration under the section `backend`. For now, we'll give your backend the name `SGE`, but you can use any name you would like.

```hocon
backend {
  providers {
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        # to be filled in
      }
    }
  }
}
```

The `actor-factory` above tells cromwell that you will be using the `config` section to tell cromwell how to submit jobs, abort jobs, etc.

You'll likely also want to change the default backend to your new backend, by setting this configuration value:

```hocon
backend.default = SGE
```

#### Specifying the runtime attributes for your HPC tasks

In the config section for your backend, you can define the different [runtime attributes](../RuntimeAttributes) that your HPC tasks will support. Any runtime attribute configured here will be read from the WDL tasks, and then passed into the command line used to submit jobs to the HPC cluster.

All runtime attributes must be defined in a single multi-line block. The syntax of this block is the same as defining the inputs for a WDL task.

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  Int cpu = 1
  Float? memory_gb
  String? sge_queue
  String? sge_project
  """
}
```

In the example above, we have defined four different WDL variables defined, `cpu`, `memory_gb`, `sge_queue`, and `sge_project`. Below you will find more information on `cpu` and `memory`, and the ability to add custom runtime attributes like the `sge_queue` and `sge_project`.

**cpu**

When you declare a runtime attribute with the name `cpu`, it must be an `Int`. This integer will validated to always be `>= 1`.

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  Int cpu = 1
  # ...
  """
}
```

**memory**

When running a workflow, the memory runtime attribute in the task will specify the units of memory. For example, this jobs specifies that it only needs 512 megabytes of memory when running.

```wdl
task hello {
  command { echo hello }
  runtime { memory: "512 MB" }
}
```

However, it's possible that when submitting jobs to your HPC cluster you want to specify the units in gigabytes.

To specify the memory units that the submit command should use, append the units to the memory runtime attribute. For example:

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  Float? memory_gb
  # ...
  """
}
```

Now, no matter what unit of memory is used within the task, the value will be converted into gigabytes before it is passed to your submit command.

**custom attributes**

You can also declare other runtime attributes that a WDL task may use. For example, suppose you would like to allow the WDL to specify an sge queue in a task, like:

```wdl
task hello {
  command { echo hello }
  runtime { sqe_queue: "short" }
}
```

You declare your runtime attribute in your config by adding a new input:

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  String? sge_queue
  # ...
  """
}
```

In this case, we've stated that the `sge_queue` is optional. This allows us to reuse WDLs from other pipeline authors who may not have set an `sge_queue`.

Alternatively, you can also set a default for the declared runtime attributes.

```hocon
backend.providers.SGE.config {
  runtime-attributes = """
  String sge_queue = "short"
  # ...
  """
}
```

#### How Cromwell should start an HPC job

When Cromwell runs a task, it will fill in a template for the job using the declared runtime attributes. This specific template will vary depending on the requirements of your HPC cluster. For example, say you normally submit jobs to SGE using:

```bash
qsub -terse -V -b y -N my_job_name \
  -wd /path/to/working_directory \
  -o /path/to/stdout.qsub \
  -e /path/to/stderr.qsub \
  -pe smp 1 -l mem_free=0.5g -q short \
  /usr/bin/env bash myScript.bash
```

For this particular SGE cluster, the above sets the working directory, stdout and stderr paths, the number of cpus to 1, the memory to half a gigabyte, and runs on the short queue.

Converting this into a template using our runtime attributes requires defining `submit` as one would a WDL task `command`:

```hocon
backend.providers.SGE.config {
  submit = """
  qsub \
  -terse \
  -V \
  -b y \
  -N ${job_name} \
  -wd ${cwd} \
  -o ${out}.qsub \
  -e ${err}.qsub \
  -pe smp ${cpu} \
  ${"-l mem_free=" + memory_gb + "g"} \
  ${"-q " + sge_queue} \
  ${"-P " + sge_project} \
  /usr/bin/env bash ${script}
  """
}
```

When the job finishes submitting, Cromwell will need to retrieve the job id, so that it can abort the job if necessary. This job should be written to the stdout after submission, where Cromwell will then read the job id. Because the job id may be surrounded by other text, a custom regular expression should capture the actual job id. Because the submit above uses `-terse`, the job id will be the entire contents of the stdout, but should be all digits:

```hocon
backend.providers.SGE.config {
  job-id-regex = "(\\d+)"
}
```

#### How Cromwell should abort an HPC job

When aborting an HPC job, Cromwell will run a command confifured under the key `kill`, passing in the WDL variable `job_id`:

```hocon
backend.providers.SGE.config {
  kill = "qdel ${job_id}"
}
```

#### How Cromwell checks if an HPC job is alive

Whenever Cromwell restarts it checks to see if a job has completed by searching for return code in a file called `rc`. If this file isn't available, in this case Cromwell runs an extra check to make sure the job is still alive. You can configure the command used for this check via:

```hocon
backend.providers.SGE.config {
  check-alive = "qstat -j ${job_id}"
}
```

#### Other backend settings

On some systems, the administrators may limit the number of HPC jobs a user may run at a time. To configure this limit, you can use the value `concurrent-job-limit` to limit the number of jobs.

```hocon
backend.providers.SGE.config {
  concurrent-job-limit = 100
}
```

#### Putting the config section all together

With the above sections, we can combine them all together to create a completly working HPC backend.

```hocon
backend {
  default = SGE
  
  providers {
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        concurrent-job-limit = 100
    
        runtime-attributes = """
        Int cpu = 1
        Float? memory_gb
        String? sge_queue
        String? sge_project
        """
    
        submit = """
        qsub \
        -terse \
        -V \
        -b y \
        -N ${job_name} \
        -wd ${cwd} \
        -o ${out} \
        -e ${err} \
        -pe smp ${cpu} \
        ${"-l mem_free=" + memory_gb + "g"} \
        ${"-q " + sge_queue} \
        ${"-P " + sge_project} \
        /usr/bin/env bash ${script}
        """

        job-id-regex = "(\\d+)"

        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
      }
    }
  }
}
```

Running Cromwell with this in our configuration file will now submit jobs to SGE!

### Next steps

You might find the following tutorials interesting to tackle next:

* [Persisting Data Between Restarts](PersistentServer)
* [Server Mode](ServerMode.md)
