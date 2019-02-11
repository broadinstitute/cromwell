## Containers

Containers are technologies that describe an environment's software requirements. Broadly, there are two primary types of containers:

* Images - Software requirements that describe a sandbox environment that you.
* Runtime - Specifies a series of processes that run inside a container image.

Containers aim to be content-agnostic, infrastructure-agnostic, designed for automation and promote reproducibility.

There are organisations such as the [Open Container Initiative](https://www.opencontainers.org) that are aiming to standardise the interface that container technologies can implement.

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

### Let's get started

Container's are specified on a per-task level, this can be achieved in WDL by specifying a [`docker`](https://software.broadinstitute.org/wdl/documentation/spec#docker) tag in the `runtime` section. For example, specifying that the following WDL script should use the container `ubuntu:latest` can be achieved by:

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

<!-- Similar in CWL, you can specifiy a [`DockerRequirement`](https://www.commonwl.org/v1.0/CommandLineTool.html#DockerRequirement) inside the requirements section:

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

#### Singularity

#### udocker

### Enforcing container requirements

You can enforce container requirements by configuring Cromwell to only accept `submit-docker` strings.


### Notes

#### Docker Digests

_Notes about disabling docker digests_
 