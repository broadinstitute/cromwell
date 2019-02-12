## Containers

Containers are self-contained software archives, that hold everything from the tool you want to run, to libraries and runtimes it requires, and the operating system it runs on.
 
Best Practise WDL and CWL define containers for their tasks to run in, to ensure reproducibility and portability - that running the same task on a different system will run the exact same software. 

Docker images are the most common container format, but is is not advisable for certain systems to run Docker itself, and for this reason Cromwell supports a number of alternatives

### Prerequisites
This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)
* [Configuration Files](ConfigurationFiles.md)

### Goals

At the end of this tutorial, you'll become familiar with container technologies and how to configure Cromwell to use these indepdently, or with job schedulers.

We'll discuss:

* [Docker](https://www.docker.com)
* [Singularity](https://www.sylabs.io/docs/)
* [udocker](https://github.com/indigo-dc/udocker)

### Specifying Containers in your Workflow

Containers are specified on a per-task level, this can be achieved in WDL by specifying a [`docker`](https://software.broadinstitute.org/wdl/documentation/spec#docker) tag in the `runtime` section. For example, specifying that the following WDL script should use the container `ubuntu:latest` can be achieved by:

```wdl
task hello_world {
    String name = "World"
    command {
        echo 'Hello, ${name}'
    }
    output {
        File out = stdout()
    }
    runtime {
        docker: 'ubuntu:latest'
    }
}

workflow hello {
    call hello_world
}
```

Similarly in CWL, you can specify a [`DockerRequirement`](https://www.commonwl.org/v1.0/CommandLineTool.html#DockerRequirement) inside the requirements section:

```cwl
cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
    name:
        type: string
        default: "World"
        inputBinding:
          prefix: "Hello, "
outputs:
    out: stdout

requirements:
    DockerRequirement:
        dockerPull: "ubuntu:latest"
``` -->

#### Docker

Docker is a popular container technology that is natively supported by Cromwell and WDL. No extra configuration must be provided to allow docker to run (provided Docker is installed).

##### Shortfalls of Docker

Docker can allow running users to gain superuser privileges, called the [Docker daemon attack surface](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface). In HPC and multi-user environments, Docker recommends that "only trusted users should be allowed to control your Docker Daemon" [source](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface). For this reason it's worth exploring other technologies that will support the reproducibility and simplicity of running a workflow that's tagged with docker containers.

#### Docker on GCP and AWS

Docker is a required tag on some cloud providers as it is used to request or provision resources. It's currently unknown whether you could run other container technologies on these cloud providers however given the same docker image, they should give exactly the same results.

#### Singularity

<!-- Singularity introduction here -->

##### Local environments

In local environments, you must configure Cromwell to use a different `submit-docker` script that would start Singularity instead of docker. Singularity requires docker images to be prefixed with the prefix `docker://`. The submit string for Singularity is:
```bash
  singularity exec --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${script}
```

As the Singularity container does not emit a job-id, we must include the `run-in-background` tag within the the provider section in addition to the docker-submit script. As Cromwell watches for the existence of the `rc` file, the `run-in-background` option has the caveat that we require the Singularity container to successfully complete, otherwise the workflow might hang indefinitely.

Putting this together, we have an example base configuration for a local environment:
```hocon
include required(classpath("application"))

backend {
    default: singularity
    providers: {
        singularity {
            # The backend custom configuration.
            actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"

            config {
                run-in-background = true
                runtime-attributes = """
                  String? docker
                """
                submit-docker = """
                  singularity exec --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${script}
                """
            }
        }
    }
}
```

##### Job schedulers

The premise when running Singularity on job schedulers it include the singularity submit script as part of the wrapped command. On some HPCs, worker nodes do not have stable access to the internet or build access, for this reason within the submit-docker configuration we'll include a `singularity build` command that is run on the worker node. In the following example configuration for SLURM, we set a cache directory for Singularity so that once we've built the image, the worker nodes are able to access the image without any work.

The following configuration also includes the `--userns` argument, which will "run \[the\] container in a new user namespace, allowing Singularity to run completely unprivileged on recent kernels" (source: Singuarlity 3.0.3 CLI).

```hocon
include required(classpath("application"))

backend {
  default: SLURM
  providers: {
    SLURM {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
          Int runtime_minutes = 60
          Int cpus = 2
          Int requested_memory_mb_per_core = 8000
          String? queue
          String? docker
        """
        submit-docker = """
          export SINGULARITY_CACHEDIR=/path/to/singularity_cache
          # ensure singularity is loaded, potentially by: module load Singularity/3.0.1
          IMAGE=/path/to/docker_location/${docker}
          singularity build --sandbox $IMAGE docker://${docker} > /dev/null
          sbatch -J ${job_name} -D ${cwd} -o ${cwd}/execution/stdout -e ${cwd}/execution/stderr ${"-p " + queue} \
            -t ${runtime_minutes} ${"-c " + cpus} --mem-per-cpu=${requested_memory_mb_per_core} \
            --wrap "singularity exec --userns -B ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${script}"
          """
        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}
```

#### udocker

### Enforcing container requirements

You can enforce container requirements by not including the standard `submit` attribute on the provider attributes. You should only use the `submit-docker` strings, however note that some environment variables (`stdout`, `stderr`) are different between these two.

### Notes

#### Job signalling

Cromwell uses the presence of an `rc` file in the execution directory 

#### Recommendations

It's recommended against using `:latest` in favour of versioned tags or the more unique digests, as the latest tag can change through time, reducing the chance that your workflow will work or is exactly reproducible.

#### Docker Digests

_Notes about disabling docker digests_
 