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

Container's are specified on a per-task level, this can be achieved in WDL ou specify a container in WDL as a Docker tag in the runtime section of a task

#### Docker

Docker is a popular container technology that is natively supported by Cromwell and WDL. No extra configuration must be provided to allow docker to run (provided Docker is installed).

##### Shortfalls of Docker

#### Singularity

#### udocker

### Enforcing container requirements

You can enforce container requirements by configuring Cromwell to only accept `submit-docker` strings.


### Notes

#### Docker Digests

_Notes about disabling docker digests_
 