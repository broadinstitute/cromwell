**HTCondor Backend**

Allows to execute jobs using HTCondor which is a specialized workload management system for compute-intensive jobs created by the Center for High Throughput Computing in the Department of Computer Sciences at the University of Wisconsin-Madison (UW-Madison).

The backend is specified via the actor factory `ConfigBackendLifecycleActorFactory`:

```
backend {
  providers {
    HtCondor {
      config {
        actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
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

The job is launched with a configurable script command such as:

```
chmod 755 ${script}
cat > ${cwd}/execution/submitFile <<EOF
Iwd=${cwd}/execution
requirements=${nativeSpecs}
leave_in_queue=true
request_memory=${memory_mb}
request_disk=${disk_kb}
error=${err}
output=${out}
log_xml=true
request_cpus=${cpu}
executable=${script}
log=${cwd}/execution/execution.log
queue
EOF
condor_submit ${cwd}/execution/submitFile
```

The HtCondor backend gets the job ID from parsing the `submit.stdout` text file.

Since the `script.sh` ends with `echo $? > rc`, the backend will wait for the existence of this file, parse out the return code and determine success or failure and then subsequently post-process.

The command used to submit the job is specified under the configuration key `backend.providers.HtCondor.config.submit`. It uses the same syntax as a command in WDL, and will be provided the variables:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `cwd` - The path where the script should be run.
* `out` - The path to the stdout.
* `err` - The path to the stderr.
* `job_name` - A unique name for the job.

This backend also supports docker as optional feature. Configuration key `backend.providers.HtCondor.config.submit-docker` is specified for this end. When the WDL contains a docker runtime attribute, this command will be provided with two additional variables:

* `docker` - The docker image name.
* `docker_cwd` - The path where `cwd` should be mounted within the docker container.
* `docker_script` - The path of the `script` within the docker container.
* `docker_out` - The path of the `out` within the docker container.
* `docker_err` - The path of the `err` within the docker container.

```
chmod 755 ${script}
cat > ${cwd}/execution/dockerScript <<EOF
#!/bin/bash
docker run --rm -i -v ${cwd}:${docker_cwd} ${docker} /bin/bash ${docker_script}
EOF
chmod 755 ${cwd}/execution/dockerScript
cat > ${cwd}/execution/submitFile <<EOF
Iwd=${cwd}/execution
requirements=${nativeSpecs}
leave_in_queue=true
request_memory=${memory_mb}
request_disk=${disk_kb}
error=${cwd}/execution/stderr
output=${cwd}/execution/stdout
log_xml=true
request_cpus=${cpu}
executable=${cwd}/execution/dockerScript
log=${cwd}/execution/execution.log
queue
EOF
condor_submit ${cwd}/execution/submitFile
```

This backend support additional runtime attributes that are specified in the configuration key `backend.providers.HtCondor.config.runtime-attributes`. It uses the same syntax as specifying runtime attributes in a task in WDL.

There are five special runtime attribute configurations, `cpu`, `memory_mb`, `disk_kb`, `nativeSpecs`, `docker`.
Optional values are defined with the prefix `?` attached to the type.

```
backend {
  providers {
    HtCondor {
      config {
        # ... other configuration
	    runtime-attributes = """
	       Int cpu = 1
	       Float memory_mb = 512.0
	       Float disk_kb = 256000.0
	       String? nativeSpecs
	       String? docker
	    """
      }
    }
  }
}
```

**Native Specifications**

The use of runtime attribute 'nativeSpecs' allows to the user to attach custom HtCondor configuration to tasks.
An example of this is when there is a need to work with 'requirements' or 'rank' configuration.

```
"runtimeAttributes": {
    cpu = 2
    memory = "1GB"
    disk = "1GB"
    nativeSpecs: "TARGET.Arch == \"INTEL\" && TARGET.Memory >= 64"
}
```

`nativeSpecs` attribute needs to be specified as String.

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
