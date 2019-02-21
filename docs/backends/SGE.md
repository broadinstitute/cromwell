**Sun GridEngine Backend**

The GridEngine and similar backends use programs such as `qsub` to launch a job and will poll the filesystem to determine if a job is completed.

The backend is specified via the actor factory `ConfigBackendLifecycleActorFactory`:

```
backend {
  providers {
    SGE {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        # ... other configuration
      }
    }
  }
}
```

This backend makes the same assumption about the filesystem that the local backend does: the Cromwell process and the jobs both have read/write access to the CWD of the job.

The CWD will contain a `script.sh` file which will contain the same contents as the Local backend:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

The job is launched with a configurable command such as:

```bash
qsub \
    -terse \
    -V \
    -b n \
    -N ${job_name} \
    -wd ${cwd} \
    -o ${out}.qsub \
    -e ${err}.qsub \
    -pe smp ${cpu} \
    ${"-l m_mem_free=" + memory_gb + "gb"} \
    ${"-q " + sge_queue} \
    ${"-P " + sge_project} \
    ${script}
```

The SGE backend gets the job ID from parsing the `submit.stdout` text file.

Since the `script.sh` ends with `echo $? > rc`, the backend will wait for the existence of this file, parse out the return code and determine success or failure and then subsequently post-process.

The command used to submit the job is specified under the configuration key `backend.providers.SGE.config.submit`. It uses the same syntax as a command in WDL, and will be provided the variables:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `cwd` - The path where the script should be run.
* `out` - The path to the stdout.
* `err` - The path to the stderr.
* `job_name` - A unique name for the job.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        submit = """
        qsub \
            -terse \
            -V \
            -b n \
            -N ${job_name} \
            -wd ${cwd} \
            -o ${out}.qsub \
            -e ${err}.qsub \
            ${script}
        """
      }
    }
  }
}
```

If the backend supports docker, the optional configuration keys `backend.providers.<backend>.config.submit-docker`
and  `backend.providers.<backend>.config.kill-docker` may be specified. When the WDL contains a docker runtime
attribute, this command will be provided three additional variables:

* `docker` - The docker image name.
* `docker_cwd` - The path where `cwd` should be mounted within the docker container.
* `docker_cid` - The host path to which the [container ID file](https://docs.docker.com/engine/reference/run/#pid-equivalent) should be written.
* `docker_script` - The path of the `script` inside the docker container.
* `docker_out` - The path of the `out` inside the docker container.
* `docker_err` - The path of the `err` inside the docker container.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        submit-docker = """
        qsub \
            -terse \
            -V \
            -b n \
            -N ${job_name} \
            -wd ${cwd} \
            -o ${out}.qsub \
            -e ${err}.qsub \
            -l docker,docker_images="${docker}"
            -xdv ${cwd}:${docker_cwd}
            ${script}
        """
      }
    }
  }
}
```

If the backend would like to support additional runtime attributes they may be specified in the configuration key `backend.providers.<backend>.config.runtime-attributes`. It uses the same syntax as specifying runtime attributes in a task in WDL.

There are two special runtime attribute configurations, `cpu`, and `memory_<unit>`.

When the runtime attribute configuration `Int cpu` is specified, it is always validated as a positive integer.

When the runtime attribute configuration `Int memory_<unit>` or `Float memory_<unit>` is specified, it is provided to submit by the runtime attribute in WDL `memory`.

For example, if the backend specifies the configuration for `backend.providers.<backend>.config.runtime-attributes` as:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = "Float memory_mb"
      }
    }
  }
}
```

And the WDL specifies a task with:

```
task hello_gigabyte {
  command { echo "hello world" }
  runtime { memory: "1 GB" }
}
```

Then for this call, the backend will be provided an additional variable `memory_mb` set to `1000.0`.

Other runtime attributes may be defined by specifying them in under the runtime attributes configuration.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float memory_mb
        String sge_project
        """
      }
    }
  }
}
```

These variables will then be passed from the WDL into the submit configuration. If one would like to have a default value, just like in WDL, the configuration may specify that the value have a default. The default must match the defined type or an error will be produced.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float memory_mb = 2.0
        String sge_project = "default"
        """
      }
    }
  }
}
```

Optional values may also be used by appending `?` to the type:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float? memory_mb
        String? sge_project
        """
      }
    }
  }
}
```

The value will be passed to the submit configuration if provided, and omitted otherwise.

There are also configuration values related to how jobs are rechecked on startup and aborted.

The option is `backend.providers.<backend>.config.run-in-background`. When `true` the backend runs the submit configuration and records the unix process id (PID). To abort the job, the PID is stopped with the unix command `kill`. Upon a cromwell restart, the PID is checked via the unix command `ps` to see if it is still alive, before cromwell goes back to polling for the `rc` file.

When `backend.providers.<backend>.config.run-in-background` is `false`, the default, the backend must specify how read the job identifier from the stdout of the submit, how to kill the job, and how to check if the job is still running during a cromwell restart. These three configuration values are `job-id-regex`, `kill`, and `check-alive`, respectively:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        job-id-regex = "(\\d+)"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        """
      }
    }
  }
}
```

The `job-id-regex` should contain one capture group while matching against the whole line or stdout file. The `check-alive` should return zero if the job is still alive.

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
