## Containers

Containers are encapsulated environments that include an operating system, libraries, and software. For example, if you have a host machine running Centos, you can run an isolated container with Ubuntu 18.04. At a high level, it's useful to think of a container as a program or binary.
 
To promote reproducibility and portability, it's considered best practice to define containers for a WDL and CWL task to run in - this ensures that running the same task on a different system will run the exact same software. 

Docker images are the most common container format, but it is not advisable for certain systems to run Docker itself, and for this reason Cromwell can be configured to support a number of alternatives.

* [Prerequisites](#prerequisites)
* [Goals](#goals)
* [Specifying Containers in your Workflow](#specifying-containers-in-your-workflow)
* [Docker](#docker)
    * [Docker on a Local Backend](#docker-on-a-local-backend)
    * [Docker on Cloud](#docker-on-cloud)
    * [Docker on HPC](#docker-on-hpc)
* [Singularity](#singularity)
    * [Installation](#installation)
    * [Configuring Cromwell for Singularity](#configuring-cromwell-for-singularity)
        * [Local environments](#local-environments)
        * [Job schedulers](#job-schedulers)
    * [Without Setuid](#without-setuid)
    * [Singularity Cache](#singularity-cache)
* [udocker](#udocker)
    * [Installation](#installation-1)
    * [Configuration](#configuration)
    * [Caching](#caching)
* [Configuration in Detail](#configuration-in-detail)
    * [Enforcing container requirements](#enforcing-container-requirements)
    * [Docker Digests](#docker-digests)
    * [Docker Root](#docker-root)
    * [Docker Config Block](#docker-config-block)
* [Best Practices](#best-practices)
    * [Image Versions](#image-versions)
* [Notes](#notes)
    * [How does Cromwell know when a job or container has completed?](#how-does-cromwell-know-when-a-job-or-container-has-completed)
    * [Cromwell: Run-in-background](#cromwell-run-in-background)
* [Next Steps](#next-steps)

### Prerequisites
This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)
* [Configuration Files](ConfigurationFiles.md)
* Recommended: [Getting started on HPC clusters](HPCIntro.md)

### Goals

At the end of this tutorial, you'll become familiar with container technologies and how to configure Cromwell to use these independently, or with job schedulers.


### Specifying Containers in your Workflow

Containers are specified on a per-task level, this can be achieved in WDL by specifying a [`docker`](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#docker) tag in the `runtime` section. For example, the following script should run in the `ubuntu:latest` container:

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
``` 
___
### Docker

[Docker](https://www.docker.com) is a popular container technology that is natively supported by Cromwell and WDL.


#### Docker on a Local Backend

On a single machine (laptop or server), no extra configuration is needed to allow docker to run, provided Docker is installed.

You can install Docker for Linux, Mac or Windows from [Docker Hub](https://hub.docker.com/search/?type=edition&offering=community)

#### Docker on Cloud

It is strongly advised that you provide a Docker image to tasks that will run on Cloud backends, and in fact most Cloud providers require it.

It might be possible to use an alternative container engine, but this is not recommended if Docker is supported.

#### Docker on HPC

Docker can allow running users to gain superuser privileges, called the [Docker daemon attack surface](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface). In HPC and multi-user environments, Docker recommends that "only trusted users should be allowed to control your Docker Daemon".

For this reason, this tutorial will also explore other technologies that support the reproducibility and simplicity of running a workflow that use docker containers; Singularity and udocker.

___

### Singularity

Singularity is a container technology designed for use on HPC systems in particular, while ensuring an appropriate level of security that Docker cannot provide.

#### Installation
Before you can configure Cromwell on your HPC system, you will have to install Singularity, which is documented [here](https://www.sylabs.io/guides/3.0/admin-guide/admin_quickstart.html#installation).
In order to gain access to the full set of features in Singularity, it is strongly recommended that Singularity is installed by root, with the `setuid` bit enabled, as is ([documented here](https://www.sylabs.io/guides/2.6/admin-guide/security.html#how-does-singularity-do-it)).
This likely means that you will have to ask your sysadmin to install it for you.
Because `singularity` ideally needs `setuid`, your admins may have some qualms about giving Singularity this privilege.
If that is the case, you might consider forwarding [this letter](https://www.sylabs.io/guides/3.0/user-guide/installation.html#singularity-on-a-shared-resource) to your admins.

If you are not able to get Singularity installed with these privileges, you can attempt a user install.
If this is the case, you will have to alter your Cromwell configuration to work in "sandbox" mode, which is explained in [this part](#without-setuid) of the documentation. 

#### Configuring Cromwell for Singularity

Once Singularity is installed, you'll need to modify the `config` block inside `backend.providers` in your Cromwell configuration. In particular, this block contains a key called `submit-docker`, which will contain a script that is run whenever a job needs to run that uses a Docker image. If the job does not specify a Docker image, the regular `submit` block will be used.

As the configuration will require more knowledge about your execution environment, see the local and job scheduler sections below for example configurations.

##### Local environments

On local backends, you have to configure Cromwell to use a different 
`submit-docker` script that would start Singularity instead of docker. 
Singularity requires docker images to be prefixed with the prefix `docker://`.

Using containers isolates the filesystem that the script is allowed to interact with, for that reason we'll bind in the current working directory as `${docker_cwd}`, and we'll use the container-specific script path `${docker_script}`. 

An example submit script for Singularity is:
```bash
singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}
```

As the `Singularity exec` command does not emit a job-id, we must include the `run-in-background` tag within the the provider section in addition to the docker-submit script. As Cromwell watches for the existence of the `rc` file, the `run-in-background` option has the caveat that we require the Singularity container to successfully complete, otherwise the workflow might hang indefinitely.

To ensure reproducibility and an isolated environment inside the container, 
`--containall` is an **important** function. By default, Singularity will mount
the user's home directory and import the user's environment as well as some 
other things that make Singularity easier to use in an interactive shell. 
Unfortunately settings in the home directory and the user's environment may 
affect the outcome of the tools that are used. This means different users may
get different results. Therefore, to ensure reproducibility while using 
Singularity, the `--containall` flag should be used. This will make sure the 
environment is cleaned and the HOME directory is not mounted.

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
                  singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}
                """
            }
        }
    }
}
```

##### Job schedulers

To run Singularity on a job scheduler, the singularity command needs to be passed to the scheduler as a wrapped command.

For example, in SLURM, we can use the normal SLURM configuration as explained in the [SLURM documentation](../backends/SLURM), however we'll add a `submit-docker` block to execute when a task is tagged with a docker container. 

When constructing this block, there are a few things to keep in mind:
- Make sure Singularity is loaded (and in PATH). If `module` is installed for 
  example you can call `module load Singularity`. If the cluster admin has made
  a Singularity module available. Alternatively you can alter the `PATH` 
  variable directly or simply use `/path/to/singularity` 
  directly in the config.
- We should treat worker nodes as if they do not have stable access to the 
  internet or build access, so we will pull the container before the task is 
  submit to the cluster.
- It's a good idea to use a Singularity cache so that same images should only
  have to be pulled once. Make sure you set the `SINGULARITY_CACHEDIR` 
  environment variable to a location on the filesystem that is reachable by the
  worker nodes!
- If we are using a cache we need to ensure that submit processes started by
  Cromwell do not pull to the same cache at the same time. This may corrupt the
  cache. We can prevent this by implementing a filelock with `flock` and 
  pulling the image before the job is submitted. The flock and pull command 
  needs to be placed *before* the submit command so all pull commands are 
  executed on the same node. This is necessary for the filelock to work.
- As mentioned above the `--containall` flag is **important** for 
  reproducibility.

```
submit-docker = """
    # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
    # based on the users home.
    if [ -z $SINGULARITY_CACHEDIR ]; 
        then CACHE_DIR=$HOME/.singularity/cache
        else CACHE_DIR=$SINGULARITY_CACHEDIR
    fi
    # Make sure cache dir exists so lock file can be created by flock
    mkdir -p $CACHE_DIR  
    LOCK_FILE=$CACHE_DIR/singularity_pull_flock
    # Create an exclusive filelock with flock. --verbose is useful for 
    # for debugging, as is the echo command. These show up in `stdout.submit`.
    flock --verbose --exclusive --timeout 900 $LOCK_FILE \
    singularity exec --containall docker://${docker} \
    echo "successfully pulled ${docker}!"

    # Submit the script to SLURM
    sbatch \
      [...]
      --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
  """
```

Putting this all together, a complete SLURM + Singularity config might look like this: 

```
backend {
  default = slurm

  providers {
    slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"                                                                                     
      config {
        runtime-attributes = """
        Int runtime_minutes = 600
        Int cpus = 2
        Int requested_memory_mb_per_core = 8000
        String? docker
        """

        submit = """
            sbatch \
              --wait \
              -J ${job_name} \
              -D ${cwd} \
              -o ${out} \
              -e ${err} \
              -t ${runtime_minutes} \
              ${"-c " + cpus} \
              --mem-per-cpu=${requested_memory_mb_per_core} \
              --wrap "/bin/bash ${script}"
        """

        submit-docker = """
            # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
            # based on the users home.
            if [ -z $SINGULARITY_CACHEDIR ]; 
                then CACHE_DIR=$HOME/.singularity/cache
                else CACHE_DIR=$SINGULARITY_CACHEDIR
            fi
            # Make sure cache dir exists so lock file can be created by flock
            mkdir -p $CACHE_DIR  
            LOCK_FILE=$CACHE_DIR/singularity_pull_flock
            # Create an exclusive filelock with flock. --verbose is useful for 
            # for debugging, as is the echo command. These show up in `stdout.submit`.
            flock --verbose --exclusive --timeout 900 $LOCK_FILE \
            singularity exec --containall docker://${docker} \
            echo "successfully pulled ${docker}!"

            # Submit the script to SLURM
            sbatch \
              --wait \
              -J ${job_name} \
              -D ${cwd} \
              -o ${cwd}/execution/stdout \
              -e ${cwd}/execution/stderr \
              -t ${runtime_minutes} \
              ${"-c " + cpus} \
              --mem-per-cpu=${requested_memory_mb_per_core} \
              --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
        """

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}

```

#### Without Setuid
In addition, if you or your sysadmins were not able to give `setuid` permissions to `singularity`, you'll have to modify the config further to ensure the use of sandbox images:

```
submit-docker = """
    [...]

    # Build the Docker image into a singularity image
    # We don't add the .sif file extension because sandbox images are directories, not files
    DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
    IMAGE=${cwd}/$DOCKER_NAME
    singularity build --sandbox $IMAGE docker://${docker}

    # Now submit the job
    # Note the use of --userns here
    sbatch \
      [...]
      --wrap "singularity exec --userns --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}"
"""
```

#### Singularity Cache
By default, Singularity will cache the Docker images you pull in `~/.singularity`, your home directory.

However, if you are sharing your Docker images with other users or have limited space in your user directory, you can redirect this caching location by exporting the `SINGULARITY_CACHEDIR` variable in your `.bashrc` or at the start of the `submit-docker` block.
```
export SINGULARITY_CACHEDIR=/path/to/shared/cache
```

For further information on the Singularity Cache, refer to the [Singularity 2 caching documentation](https://www.sylabs.io/guides/2.6/user-guide/build_environment.html#cache-folders) (this hasn't yet been updated for Singularity 3).

___

### udocker

[udocker](https://github.com/indigo-dc/udocker) is a tool designed to "execute simple docker containers in user space without requiring root privileges".

In essence, udocker provides a command line interface that mimics `docker`, and implements the commands using one of four different container backends:

* PRoot
* Fakechroot
* runC
* Singularity

#### Installation
udocker can be installed without any kind of root permissions. Refer to udocker's installation documentation [here](https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md) for more information.

#### Configuration

(As of [2019-02-18](https://github.com/indigo-dc/udocker/issues/112)) udocker does not support looking up docker container by digests, hence you'll have to make ensure `hash-lookup` is disabled. Refer to [this section](#docker-digests) for more detail.

To configure `udocker` to work in a local environment, you must tag the provider's configuration to `run-in-background` and update the `submit-docker` to use udocker:
```
run-in-background = true
submit-docker = """
    udocker run -v ${cwd}:${docker_cwd} ${docker} ${job_shell} ${docker_script}
"""
```

With a job queue like SLURM, you just need to wrap this script in an `sbatch` submission like we did with Singularity:

```
submit-docker = """
    # Pull the image using the head node, in case our workers don't have network access
    udocker pull ${docker}
    
    sbatch \
      -J ${job_name} \
      -D ${cwd} \
      -o ${cwd}/execution/stdout \
      -e ${cwd}/execution/stderr \
      -t ${runtime_minutes} \
      ${"-c " + cpus} \
      --mem-per-cpu=${requested_memory_mb_per_core} \
      --wrap "udocker run -v ${cwd}:${docker_cwd} ${docker} ${job_shell} ${docker_script}"
"""
```

#### Caching
udocker caches images in a single directory, which defaults to [`~/.udocker`](https://github.com/indigo-dc/udocker/blob/master/udocker.py#L137), meaning that caching is done on a per-user basis. 
However, like Singularity, if you want to share a cache with other users in your project,you you can override the location of the udocker cache directory either using:
* A config file [described here](https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md#9-configuration), containing a line such as `topdir = "/path/to/cache"`.
* Using the environment variable `$UDOCKER_DIR`

___

### Configuration in Detail
The behaviour of Cromwell with containers can be modified using a few other options.

#### Enforcing container requirements
You can enforce the use of a container by not including the `submit` block in the provider section.

However note that some interpolated variables (`${stdout}`, `${stderr}`) are different between these two blocks.

#### Docker Digests

Each Docker repository has a number of tags that can be used to refer to the latest image of a particular type. 
For instance, when you run a normal Docker image with `docker run image`, it will actually run `image:latest`, the `latest` tag of that image.

However, by default Cromwell requests and runs images using their `sha` hash, rather than using tags.
This strategy is actually preferable, because it ensures every execution of the task or workflow will use the exact same version of the image, but some engines such as `udocker` don't support this feature.

If you are using `udocker` or want to disable the use of hash-based image references, you can set the following config option:
```
docker.hash-lookup.enabled = false
```

Nb: By disabling hash-lookup, call caching will not work for any container using a floating tag.

#### Docker Root
If you want to change the root directory inside your containers, where the task places input and output files, you can edit the following option:

```
backend {
  providers {
    LocalExample {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
      
        # Root directory where Cromwell writes job results in the container. This value
        # can be used to specify where the execution folder is mounted in the container.
        # it is used for the construction of the docker_cwd string in the submit-docker
        # value above.
        dockerRoot = "/cromwell-executions"
      }
    }
  }
}
```


#### Docker Config Block
Further docker configuration options available to be put into your config file are as follows. 
For the latest list of parameters, refer to the [example configuration file][cromwell-examples-conf],
and [specific backend provider examples][cromwell-examples-folder].

```
docker {
  hash-lookup {
    # Set this to match your available quota against the Google Container Engine API
    #gcr-api-queries-per-100-seconds = 1000

    # Time in minutes before an entry expires from the docker hashes cache and needs to be fetched again
    #cache-entry-ttl = "20 minutes"

    # Maximum number of elements to be kept in the cache. If the limit is reached, old elements will be removed from the cache
    #cache-size = 200

    # How should docker hashes be looked up. Possible values are "local" and "remote"
    # "local": Lookup hashes on the local docker daemon using the cli
    # "remote": Lookup hashes on docker hub, gcr, gar, quay
    #method = "remote"
  }
}
```


### Best Practices

#### Image Versions

When choosing the image version for your pipeline stages, it is highly recommended that you use a hash rather than a tag, for the sake of reproducibility
For example, in WDL, you could do this:
```wdl
runtime {
    docker: 'ubuntu:latest'
}
```

But what you should do is this:
```wdl
runtime {
    docker: 'ubuntu@sha256:7a47ccc3bbe8a451b500d2b53104868b46d60ee8f5b35a24b41a86077c650210'
}
```

You can find the `sha256` of an image using `docker images --digests`
 
 
### Notes

#### How does Cromwell know when a job or container has completed?
Cromwell uses the presence of the `rc` (returncode) file to determine whether a task has succeeded or failed. This `rc` file is generated as part of the `script` within the execution directory, where the script is assembled at runtime. This is important as if the script executes successfully but the container doesn't terminate, Cromwell will continue the execution of the workflow and the container will persist hogging system resources.

Within the configurations above:
- `singularity`: The exec mode does not run a container on the background

#### Cromwell: Run-in-background

By enabling Cromwell's run-in-background mode, you remove the necessity for the `kill`, `check-alive` and `job-id-regex` blocks, which disables some safety checks when running workflows:

- If there is an error starting the container or executing the script, Cromwell may not recognise this error and hang. For example, this may occur if the container attempts to exceed its allocated resources (runs out of memory); the container daemon may terminate the container without completing the script.
- If you abort the workflow (by attempting to close Cromwell or issuing an abort command), Cromwell does not have a reference to the container execution and will not be able to terminate the container.

This is only necessary in local environments where there is no job manager to control this, however if your container technology can emit an identifier to stdout, then you are able to remove the run-in-background flag. 

### Next Steps

Congratulations for improving the reproducibility of your workflows! You might find the following cloud-based tutorials interesting to test your workflows (and ensure the same results) in a completely different environment:

- [Getting started with AWS Batch](AwsBatch101.md)
- [Getting started on Google Pipelines API](PipelinesApi101.md)

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends
