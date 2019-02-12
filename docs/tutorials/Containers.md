## Containers

Containers are self-contained software archives, that hold everything from the tool you want to run, to libraries and runtimes it requires, and the operating system it runs on.
 
Best Practise WDL and CWL define containers for their tasks to run in, to ensure reproducibility and portability - that running the same task on a different system will run the exact same software. 

Docker images are the most common container format, but it is not advisable for certain systems to run Docker itself, and for this reason Cromwell supports a number of alternatives

### Prerequisites
This tutorial page relies on completing the previous tutorials:

* [Five Minute Introduction](FiveMinuteIntro.md)
* [Configuration Files](ConfigurationFiles.md)
* Recommended: [Getting started on HPC clusters](HPCIntro.md)

### Goals

At the end of this tutorial, you'll become familiar with container technologies and how to configure Cromwell to use these independently, or with job schedulers.

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
``` 

### Docker

[Docker](https://www.docker.com) is a popular container technology that is natively supported by Cromwell and WDL.


#### Docker on a Local Backend

On a single machine (laptop or server), no extra configuration is needed to allow docker to run, provided Docker is installed.

You can install Docker for Linux, Mac or Windows from [Docker Hub](https://hub.docker.com/search/?type=edition&offering=community)

#### Docker on Cloud

It is strongly advised that you use Docker images on Cloud backends, because most cloud backends require Docker images to run.

It might be possible to use an alternative container engine, but this is not recommended if Docker is available.

#### Docker on HPC

Docker can allow running users to gain superuser privileges, called the [Docker daemon attack surface](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface). In HPC and multi-user environments, Docker recommends that "only trusted users should be allowed to control your Docker Daemon".

For this reason it's worth exploring other technologies that will support the reproducibility and simplicity of running a workflow that uses docker containers; Singularity or uDocker.

### Singularity

Singularity is a container engine designed for use on HPC systems in particular, while ensuring an appropriate level of security that Docker cannot provide.

#### Installation
Before you can configure Cromwell on your HPC system, you (or your sysadmin) will have to install Singularity, which is documented [here](https://www.sylabs.io/guides/3.0/admin-guide/admin_quickstart.html#installation).

When installing Singularity, your admins have to make a decision: do they grant `setuid` to the Singularity binary or not?
[As is documented here](https://www.sylabs.io/guides/2.6/admin-guide/security.html#how-does-singularity-do-it), Singularity can run either using "sandbox" containers, or using proper image files.
However, unless `setuid` is set, you are restricted to using sandbox images, which come with a number of downsides.
If you can't convince Singularity can be trusted with `setuid`, don't fear, Cromwell will still work.
However, this is not ideal, so you might consider forwarding [this letter](https://www.sylabs.io/guides/3.0/user-guide/installation.html#singularity-on-a-shared-resource) to your admins.

#### Configuring Cromwell for Singularity

Once Singularity is installed, the main configuration block you'll need to modify is the `config` block inside `backend.providers` in your Cromwell configuration.
In particular, this block contains a key called `submit-docker`, which will contain a script that is run whenever a job needs to run that uses a Docker image.
If the job does not specify a Docker image, then the regular `submit` block will be used instead.

##### Local environments

On local backends, you have to configure Cromwell to use a different `submit-docker` script that would start Singularity instead of docker.
Singularity requires docker images to be prefixed with the prefix `docker://`.
An example submit script for Singularity is:
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

At its most basic level, if we want to run Singularity on a job scheduler, we need to pass the singularity command to the scheduler as a wrapped command.
If we were using SLURM, this means we can use all the normal SLURM configuration as explained in the [SLURM documentation](../backends/SLURM). 
However, our submit script will be a wrapped version of the script above:
```
submit-docker = """
  # Ensure singularity is loaded if it's installed as a module
  module load Singularity/3.0.1
  
  # Submit the script to SLURM
  sbatch \
    -J ${job_name} \
    -D ${cwd} \
    -o ${cwd}/execution/stdout \
    -e ${cwd}/execution/stderr \
    ${"-p " + queue} \
    -t ${runtime_minutes} \
    ${"-c " + cpus} \
    --mem-per-cpu=${requested_memory_mb_per_core} \
    --wrap "singularity exec --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${script}"
  """
```

On some HPC systems, worker nodes do not have stable access to the internet or build access.
If this is the case, you'll need to pull the Docker image before you actually submit the SLURM job:
```
submit-docker = """
    [...]
  
    # Build the Docker image into a singularity image, using the head node
    IMAGE=/path/to/image/${docker}.sinf
    singularity build $IMAGE docker://${docker}
  
    # Submit the script to SLURM
    sbatch \
      [...]
      --wrap "singularity exec --bind ${cwd}:${docker_cwd} ${IMAGE} ${job_shell} ${script}"
  """
```

In addition, if you or your sysadmins were not able to give `setuid` permissions to `singularity`, you'll have to modify the config further to ensure the use of sandbox images:

```
submit-docker = """
    [...]
    
    # Build the Docker image into a singularity image
    # We don't add the .sinf file extension because sandbox images are directories, not files
    IMAGE=/path/to/image/${docker}
    singularity build --sandbox $IMAGE docker://${docker}
    
    # Now submit the job
    # Note the use of --userns here
    sbatch \
      [...]
      --wrap "singularity exec --userns --bind ${cwd}:${docker_cwd} ${IMAGE} ${job_shell} ${script}"
"""
```

You might also want to cache the downloaded Docker images, by setting the `SINGULARITY_CACHEDIR` variable.
This will save unnecessary network access from downloading the same image multiple times.

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
            sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${runtime_minutes} \
            ${"-c " + cpus} \
            --mem-per-cpu=${requested_memory_mb_per_core} \
            --wrap "/bin/bash ${script}"
        """

        submit-docker = """
            # Cache the docker images we're downloading
            export SINGULARITY_CACHEDIR=/path/to/singularity/cache
            
            # Ensure singularity is loaded if it's installed as a module
            module load Singularity/version
            
            # Build the Docker image into a singularity image
            # We don't add the .sinf file extension because sandbox images are directories, not files
            IMAGE=/path/to/image/${docker}
            singularity build --sandbox $IMAGE docker://${docker}

            # Submit the script to SLURM
            sbatch \
              -J ${job_name} \
              -D ${cwd} \
              -o ${cwd}/execution/stdout \
              -e ${cwd}/execution/stderr \
              ${"-p " + queue} \
              -t ${runtime_minutes} \
              ${"-c " + cpus} \
              --mem-per-cpu=${requested_memory_mb_per_core} \
              --wrap "singularity exec --userns --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${script}"
        """

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}

```

### udocker

[udocker](https://github.com/indigo-dc/udocker) is a tool designed to "execute simple docker containers in user space without requiring root privileges".

<<<<<<< HEAD
udocker is a user-installable tool that allows the execution of docker containers in non-privileged environments. Similar to Singuarlity, it will need to be available on the worker node.

As of writing this tutorial ([udocker 1.1.3](https://github.com/indigo-dc/udocker/releases/tag/v1.1.3)), there was no support for [docker digests](https://github.com/indigo-dc/udocker/issues/112). For this reason, docker hash-lookup will need to be disabled, this will disable call caching for [floating tags](https://github.com/broadinstitute/cromwell/pull/4039#issuecomment-454829890).

##### Local

For a local configuration, it is enough to use a docker-submit string to start singularity, ie:

```HOCON
include required(classpath("application"))

docker.hash-lookup.enabled = false

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
          String? docker_user
        """
        submit-docker = """
          udocker run ${"--user " + docker_user} --rm -v ${cwd}:${docker_cwd} ${docker} ${script}
        """
      }
    }
  }
}
```

##### Job Scheduler

Similar to Singularity, the premise when running udocker on job schedulers is to include a udocker run command as part of the wrapped command. As HPCs often don't have internet access on worker nodes, some work is required to allow the head node to pull a container so that the worker nodes don't have to.

If you're familiar with `udocker`, please contribute (or complete the following) configuration:

```bash
export UDOCKER_CACHEDIR=/path/to/udocker_cache
# ensure udocker is loaded or in path
IMAGE=/path/to/docker_location/${docker}
# ask udocker to pull ${docker}, and store in UDOCKER_CACHEDIR
udocker pull ${docker}
sbatch -J ${job_name} -D ${cwd} -o ${cwd}/execution/stdout -e ${cwd}/execution/stderr ${"-p " + queue} \
  -t ${runtime_minutes} ${"-c " + cpus} --mem-per-cpu=${requested_memory_mb_per_core} \
  --wrap "udocker run ${"--user " + docker_user} --rm -v ${cwd}:${docker_cwd} ${docker} ${script}"
```

##### Limitations

udocker will not be able to perform tasks that require elevated privileges, some of these operations are provided in their [documentation](https://github.com/indigo-dc/udocker#Limitations).

### Enforcing container requirements
=======
In essence, udocker provides a command line interface that mimics `docker`, and implements the commands using one of four different container backends:
* PRoot
* Fakechroot
* runC
* Singularity
>>>>>>> 0f38ef79a... Elaborate on uDocker, advanced config, best practices

#### Installation
udocker can be installed without any kind of root permissions. Refer to the installation documentation [here](https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md).

#### Configuration

`udocker` mimics the `docker` command-line interface, so all you need to do to configure it is set the `submit-docker` configuration to something like this:
```
submit-docker = """
    udocker run --rm -v ${cwd}:${docker_cwd} ${docker} ${script}
"""
```

<<<<<<< HEAD
Cromwell uses the presence of an `rc` file in the execution directory or the contents of the `stderr` to determine a task's status. If you're job fails before starting, or during the startup of a container, it's likely that Cromwell will hang and fail to correctly report this error.
=======
With a job queue like SLURM, you just need to wrap this script in an `sbatch` submission like we did with Singularity:
>>>>>>> 0f38ef79a... Elaborate on uDocker, advanced config, best practices

```
submit-docker = """
    sbatch \
      -J ${job_name} \
      -D ${cwd} \
      -o ${cwd}/execution/stdout \
      -e ${cwd}/execution/stderr \
      ${"-p " + queue} \
      -t ${runtime_minutes} \
      ${"-c " + cpus} \
      --mem-per-cpu=${requested_memory_mb_per_core} \
      --wrap "udocker run --rm -v ${cwd}:${docker_cwd} ${docker} ${script}"
"""
```

In addition, to run `udocker`, you'll have to make sure `hash-lookup` is disabled.
Refer to [this section](#docker-digests) for more detail.

### Configuration in Detail
The behaviour of Cromwell with containers can be modified using a few other options.

#### Enforcing container requirements
You can enforce container requirements by not including the standard `submit` attribute on the provider attributes.
You should only use the `submit-docker` strings, however note that some environment variables (`stdout`, `stderr`) are different between these two.

#### Docker Digests

Each Docker repository has a number of tags that can be used to refer to the latest image of a particular type. 
For instance, when you run a normal Docker image with `docker run image`, it will actually run `image:latest`, the `latest` tag of that image.
However, by default Cromwell requests and runs images using their `sha` hash, rather than using tags.
This strategy is actually preferable, because it ensures every execution of the task or workflow will use the exact same version of the image, but some engines such as `udocker` don't support this feature.

If you are using `udocker` or want to enable the use of hash-based image references, you can set the following config option:
```
docker.hash-lookup.enabled = false
```

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
For the latest list of parameters, refer to [this example configuration](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.examples.conf)
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
    # "remote": Lookup hashes on docker hub and gcr
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
 