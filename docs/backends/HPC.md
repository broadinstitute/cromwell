Cromwell provides a generic way to configure a backend relying on most High Performance Computing (HPC) frameworks, and with access to a shared filesystem.

The two main features that are needed for this backend to be used are a way to submit a job to the compute cluster and to get its status through the command line.
You can find example configurations for a variety of those backends here:

* [SGE](SGE)
* [LSF](LSF)
* [SLURM](SLURM)
* [HTCondor](HTcondor)

## FileSystems

### Shared FileSystem
HPC backends rely on being able to access and use a shared filesystem to store workflow results.

Cromwell is configured with a root execution directory which is set in the configuration file under `backend.providers.<backend_name>.config.root`.  This is called the `cromwell_root` and it is set to `./cromwell-executions` by default.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.

When Cromwell runs a workflow, it first creates a directory `<cromwell_root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.
Any input files to a call need to be localized into the `<call_dir>/inputs` directory. There are different localization strategies that Cromwell will try until one works:

* `hard-link` - This will create a hard link to the file
* `soft-link` - Create a symbolic link to the file. This strategy is not enabled by default for tasks which specify a 
  Docker image and will be ignored.
* `copy` - Make a copy the file
* `cached-copy` An experimental feature. This copies files to a file cache in 
`<workflow_root>/cached-inputs` and then hard links them in the `<call_dir>/inputs` directory. 

`cached-copy` is intended for a shared filesystem that runs on multiple physical disks, where docker containers are used. 
Hard-links don't work between different physical disks and soft-links don't work with docker by default. Copying uses a lot of
space if a multitude of tasks use the same input. `cached-copy` copies the file only once to the physical disk containing
the `<workflow_root>` and then uses hard links for every task that needs the input file. This can save a lot of space.

The default order in `reference.conf` is `hard-link`, `soft-link`, `copy`

Shared filesystem localization is defined in the `config` section of each backend. The default stanza for the Local and HPC backends looks like this:

```
filesystems {
 local {
   localization: [
	 "hard-link", "soft-link", "copy"
   ]
 }
}
```

### Optional docker soft links

By default when Cromwell runs a local container it only mounts the workflow's execution directory. Thus any symbolic or
soft links pointing to files outside of the execution directory will resolve to paths that are not accessible within the
container.

As discussed above regarding `cache-copy`, `soft-link` is disabled by default on docker and other container
environments, and hard-links do not work across different physical disks.

However, it is possible to manually configure Cromwell to mount input paths such that soft links resolve outside and
inside containers.

```hocon
backend {
  default = "SlurmDocker"
  providers {
    SlurmDocker {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
        String? docker
        """
        # https://slurm.schedmd.com/sbatch.html
        submit-docker = """
        set -euo pipefail
        CACHE_DIR=$HOME/.singularity/cache
        mkdir -p $CACHE_DIR
        LOCK_FILE=$CACHE_DIR/singularity_pull_flock
        DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
        IMAGE=$CACHE_DIR/$DOCKER_NAME.sif
        (
          flock --verbose --exclusive --timeout 900 9 || exit 1
          if [ ! -e "$IMAGE" ]; then
            singularity build $IMAGE docker://${docker}
          fi
        ) 9>$LOCK_FILE
        sbatch \
          -J ${job_name} \
          -D ${cwd} \
          -o ${cwd}/execution/stdout \
          -e ${cwd}/execution/stderr \
          --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} --bind /mnt/one:/mnt/one:ro --bind /mnt/two:/mnt/two:ro $IMAGE ${job_shell} ${docker_script}"
        """
        # ... other configuration ...
        filesystems {
          local {
            caching.duplication-strategy = ["copy"]
            localization = ["soft-link", "copy"]
            docker.allow-soft-links: true
          }
        }
      }
    }
  }
}
```

The important parts of the example configuration above are:
* `config.filesystems.local.docker.allow-soft-links` set to `true`
* `config.submit-docker` containing `--bind /mnt/one:/mnt/one:ro --bind /mnt/two:/mnt/two:ro`

In this example the two directories `/mnt/one` and and `/mnt/two` will also be available within containers at their
original paths outside the container. So soft links pointing to paths under those directories will resolve during the
job execution. Note that if a user runs a workflow using an input file `/mnt/three/path/to/file` the job will fail
during execution as `/mnt/three` was not present inside the running container.

### Additional FileSystems

HPC backends (as well as the Local backend) can be configured to be able to interact with other type of filesystems, where the input files can be located for example.
Currently the only other filesystem supported is Google Cloud Storage (GCS). See the [Google section](GCPBatch) of the documentation for information on how to configure GCS in Cromwell.
Once you have a google authentication configured, you can simply add a `gcs` stanza in your configuration file to enable GCS:

```
backend.providers.MyHPCBackend {
  filesystems {
    gcs {
      # A reference to a potentially different auth for manipulating files via engine functions.
      auth = "application-default"
    }
  }
}
```

### Exit code timeout

If the cluster forcefully kills a job, it is unable to write its exit code anymore.
To address this the option `exit-code-timeout-seconds` can be used.
Cromwell will check the aliveness of the job with the `check-alive` script, every `exit-code-timeout-seconds` (polling).
When a job is no longer alive and another `exit-code-timeout-seconds` seconds have passed without an RC file being made, Cromwell can mark the job as failed.
If retries are enabled the job is submitted again.
This option will enable polling with the `check-alive` option, this could cause high load on whatever system `check-alive` calls.

When the option `exit-code-timeout-seconds` is **not** set cromwell will only execute the `check-alive` option after a restart of a cromwell server.

```
backend {
  providers {
    <backend name> {
      config {
        exit-code-timeout-seconds = 120
        # other config options
      }
    }
  }
}
```
