A backend represents a way to run the user's command specified in the `task` section.  Cromwell allows for backends conforming to
the Cromwell backend specification to be plugged into the Cromwell engine.  Additionally, backends are included with the
Cromwell distribution:

* Local / GridEngine / LSF / etc. - Run jobs as subprocesses or via a dispatcher.  Supports launching in Docker containers. Use `bash`, `qsub`, `bsub`, etc. to run scripts.
* Google JES - Launch jobs on Google Compute Engine through the Job Execution Service (JES).
* GA4GH TES - Launch jobs on servers that support the GA4GH Task Execution Schema (TES).
* HtCondor - Allows to execute jobs using HTCondor.
* Spark - Adds support for execution of spark jobs.

Backends are specified in the `backend` configuration block under `providers`.  Each backend has a configuration that looks like:

```hocon
backend {
  default = "Local"
  providers {
    BackendName {
      actor-factory = "FQN of BackendLifecycleActorFactory instance"
      config {
        key = "value"
        key2 = "value2"
        ...
      }
    }
  }
}
```

The structure within the `config` block will vary from one backend to another; it is the backend implementation's responsibility
to be able to interpret its configuration.

In the example below two backend types are named within the `providers` section here, so both
are available.  The default backend is specified by `backend.default` and must match the `name` of one of the
configured backends:

```hocon
backend {
  default = "Local"
  providers {
    Local {
      actor-factory = "cromwell.backend.impl.local.LocalBackendLifecycleActorFactory"
      config {
        root: "cromwell-executions"
        filesystems = {
          local {
            localization: [
              "hard-link", "soft-link", "copy"
            ]
          }
          gcs {
            # References an auth scheme defined in the 'google' stanza.
            auth = "application-default"
          }
        }
      }
    },
    JES {
      actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
      config {
        project = "my-cromwell-workflows"
        root = "gs://my-cromwell-workflows-bucket"
        maximum-polling-interval = 600
        dockerhub {
          # account = ""
          # token = ""
        }
        genomics {
          # A reference to an auth defined in the 'google' stanza at the top.  This auth is used to create
          # Pipelines and manipulate auth JSONs.
          auth = "application-default"
          endpoint-url = "https://genomics.googleapis.com/"
        }
        filesystems = {
          gcs {
            # A reference to a potentially different auth for manipulating files via engine functions.
            auth = "user-via-refresh"
          }
        }
      }
    }
  ]
}
```

## Backend Job Limits

You can limit the number of concurrent jobs for a backend by specifying the following option in the backend's config
stanza:

```
backend {
  ...
  providers {
    BackendName {
      actor-factory = ...
      config {
        concurrent-job-limit = 5
```

## Backend Filesystems

Each backend will utilize filesystems to store the directory structure of an executed workflow.  Currently, the backends and the type of filesystems that the backend use are tightly coupled.  In future versions of Cromwell, they may be more loosely coupled.

The backend/filesystem pairings are as follows:

* [Local Backend](#local-backend) and associated backends primarily use the [Shared Local Filesystem](#shared-local-filesystem).
* [JES Backend](#google-jes-backend) uses the [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem)

Note that while Local, SGE, LSF, etc. backends use the local or network filesystem for the directory structure of a workflow, they are able to localize inputs
from GCS paths if configured to use a GCS filesystem.  See [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem) for more details.

### Shared Local Filesystem

For the [local](#local-backend) and [Sun GridEngine](#sun-gridengine-backend) backends, the following is required of the underlying filesystem:

Cromwell is configured with a root execution directory which is set in the configuration file under `backend.providers.<backend_name>.config.root`.  This is called the `cromwell_root` and it is set to `./cromwell-executions` by default.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.

When Cromwell runs a workflow, it first creates a directory `<cromwell_root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.  For example, having a `stdout` and `stderr` file is common among both backends and they both write a shell script file to the `<call_dir>` as well.  See the descriptions below for details about backend-specific files that are written to these directories.

An example of a workflow output directory for a three-step workflow might look like this:

```
cromwell-executions/
└── three_step
    └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
        ├── call-cgrep
        │   ├── Users
        │   │   └── jdoe
        │   │       └── projects
        │   │           └── cromwell
        │   │               └── cromwell-executions
        │   │                   └── three_step
        │   │                       └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
        │   │                           └── call-ps
        │   │                               └── stdout
        │   ├── rc
        │   ├── script
        │   ├── stderr
        │   └── stdout
        ├── call-ps
        │   ├── rc
        │   ├── script
        │   ├── stderr
        │   └── stdout
        └── call-wc
            ├── Users
            │   └── jdoe
            │       └── projects
            │           └── cromwell
            │               └── cromwell-executions
            │                   └── three_step
            │                       └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
            │                           └── call-ps
            │                               └── stdout
            ├── rc
            ├── script
            ├── stderr
            └── stdout
```

WDL File

```wdl
task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow three_step {
  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}
```

In the above directory structure, you'll notice that the `call-cgrep` and `call-wc` sub-directories both contain a directory structure to point to the `stdout` file from the invocation of `ps`.  In these cases, that `stdout` file is a localized version of the one within `call-ps/stdout`.  By default both of those `stdout` files would be hard-links but they could also be symbolic links or copies of the file, depending on how Cromwell is configured (see below).  The directory structure is nested so deeply to avoid collisions.  For example, if either of these call invocations referenced two files called `stdout`, they'd collide if they were put into the same directory so the full directory structure is maintained.

Any input files to a call need to be localized into the `<call_dir>`.  There are a few localization strategies that Cromwell will try until one works.  Below is the default order specified in `application.conf` but this order can be overridden:

* `hard-link` - This will create a hard link (not symbolic) link to the file
* `soft-link` - Create a symbolic link to the file.  This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file

Shared filesystem localization is defined in the `config` section of each backend.  The default stanza for the Local, SGE, and associated backends looks like this:

```
filesystems {
 local {
   localization: [
	 "hard-link", "soft-link", "copy"
   ]
 }
}
```

### Google Cloud Storage Filesystem

On the JES backend the GCS (Google Cloud Storage) filesystem is used for the root of the workflow execution.
On the Local, SGE, and associated backends any GCS URI will be downloaded locally.  For the JES backend the `jes_gcs_root` [workflow option](#workflow-options) will take
precedence over the `root` specified at `backend.providers.JES.config.root` in the configuration file. Google Cloud Storage URIs are the only acceptable values for `File` inputs for
workflows using the JES backend.

## Local Backend

The local backend will simply launch a subprocess for each task invocation and wait for it to produce its rc file.

This backend creates three files in the `<call_dir>` (see previous section):

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

`<container_call_root>` would be equal to `<call_dir>` for non-Docker jobs, or it would be under `/cromwell-executions/<workflow_uuid>/call-<call_name>` if this is running in a Docker container.

When running without docker, the subprocess command that the local backend will launch is:

```
/bin/bash <script>"
```

When running with docker, the subprocess command that the local backend will launch is:

```
docker run --rm -v <cwd>:<docker_cwd> -i <docker_image> /bin/bash < <script>
```

> **NOTE**: If you are using the local backend with Docker and Docker Machine on Mac OS X, by default Cromwell can only
> run from in any path under your home directory.
>
> The `-v` flag will only work if `<cwd>` is within your home directory because VirtualBox with
> Docker Machine only exposes the home directory by default.  Any local path used in `-v` that is not within the user's
> home directory will silently be interpreted as references to paths on the VirtualBox VM.  This can manifest in
> Cromwell as tasks failing for odd reasons (like missing RC file)
>
> See https://docs.docker.com/engine/userguide/dockervolumes/ for more information on volume mounting in Docker.

## Google JES Backend

Google JES (Job Execution Service) is a Docker-as-a-service from Google.

### Configuring Google Project

You'll need the following things to get started:

* A Google Project (Manage/create projects [here](https://console.developers.google.com/project))
* A Google Cloud Storage bucket (View/create buckets in your project [here](https://console.cloud.google.com/storage/browser))

On your Google project, open up the [API Manager](https://console.developers.google.com/apis/library) and enable the following APIs:

* Google Compute Engine
* Google Cloud Storage
* Genomics API

If your project is `my-project` your bucket is `gs://my-bucket/`, then update your [Cromwell configuration file](#configuring-cromwell) as follows:

```hocon
backend {
  default = "JES"
  providers {
    JES {
      actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
      config {
        project = "my-project"
        root = "gs://my-bucket"
        genomics-api-queries-per-100-seconds = 1000
        .
        .
        .
      }
    }
  ]
}
```

If your project has API quotas other than the defaults set the `genomics-api-queries-per-100-seconds` value to be the lesser of the `Queries per 100 seconds per user` and `Queries per 100 seconds` quotas. This value will be used to help tune Cromwell's rate of interaction with JES.

### Configuring Authentication

The `google` stanza in the Cromwell configuration file defines how to authenticate to Google.  There are four different
authentication schemes that might be used:

* `application_default` - (default, recommended) Use [application default](https://developers.google.com/identity/protocols/application-default-credentials) credentials.
* `service_account` - Use a specific service account and key file (in PEM format) to authenticate.
* `user_account` - Authenticate as a user.
* `refresh_token` - Authenticate each individual workflow using a refresh token supplied in the workflow options.
* `user_service_account` - Authenticate each individual workflow using service account credentials supplied in the workflow options.

The `auths` block in the `google` stanza defines the authorization schemes within a Cromwell deployment:

```hocon
google {
  application-name = "cromwell"
  auths = [
    {
      name = "application-default"
      scheme = "application_default"
    },
    {
      name = "user-via-refresh"
      scheme = "refresh_token"
      client-id = "secret_id"
      client-secret = "secret_secret"
    },
    {
      name = "service-account"
      scheme = "service_account"
      service-account-id = "my-service-account"
      pem-file = "/path/to/file.pem"
    },
    {
      name = "user-service-account"
      scheme = "user_service_account"
    }
  ]
}
```

These authorization schemes can be referenced by name within other portions of the configuration file.  For example, both
the `genomics` and `filesystems.gcs` sections within a JES configuration block must reference an auth defined in this block.
The auth for the `genomics` section governs the interactions with JES itself, while `filesystems.gcs` governs the localization
of data into and out of GCE VMs.

#### Application Default Credentials

By default, application default credentials will be used.  There is no configuration required for application default
credentials, only `name` and `scheme` are required.

To authenticate, run the following commands from your command line (requires [gcloud](https://cloud.google.com/sdk/gcloud/)):

```
$ gcloud auth login
$ gcloud config set project my-project
```

#### Service Account

First create a new service account through the [API Credentials](https://console.developers.google.com/apis/credentials) page.  Go to **Create credentials -> Service account key**.  Then in the **Service account** dropdown select **New service account**.  Fill in a name (e.g. `my-account`), and select key type of JSON.

Creating the account will cause the JSON file to be downloaded.  The structure of this file is roughly like this (account name is `my-account`):

```
{
  "type": "service_account",
  "project_id": "my-project",
  "private_key_id": "OMITTED",
  "private_key": "-----BEGIN PRIVATE KEY-----\nBASE64 ENCODED KEY WITH \n TO REPRESENT NEWLINES\n-----END PRIVATE KEY-----\n",
  "client_email": "my-account@my-project.iam.gserviceaccount.com",
  "client_id": "22377410244549202395",
  "auth_uri": "https://accounts.google.com/o/oauth2/auth",
  "token_uri": "https://accounts.google.com/o/oauth2/token",
  "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
  "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/my-account%40my-project.iam.gserviceaccount.com"
}
```

Most importantly, the value of the `client_email` field should go into the `service-account-id` field in the configuration (see below).  The
`private_key` portion needs to be pulled into its own file (e.g. `my-key.pem`).  The `\n`s in the string need to be converted to newline characters.

While technically not part of Service Account authorization mode, one can also override the default service account that the compute VM is started with via the configuration option `JES.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`.  The service account you provide must have been granted Service Account Actor role to Cromwell's primary service account. As this only affects Google Pipelines API and not GCS, it's important that this service account, and the service account specified in `JES.config.genomics.auth` can both read/write the location specified by `JES.config.root`

#### Refresh Token

A **refresh_token** field must be specified in the [workflow options](#workflow-options) when submitting the job.  Omitting this field will cause the workflow to fail.

The refresh token is passed to JES along with the `client-id` and `client-secret` pair specified in the corresponding entry in `auths`.

#### User Service Account

A [JSON key file for the service account](####service-acocunt) must be passed in via the **user_service_account_json** field in the [workflow options](#workflow-options) when submitting the job. Omitting this field will cause the workflow to fail. The JSON should be passed as a string and will need to have no newlines and all instances of *"* and *\n* escaped. 

In the likely event that this service account does not have access to Cromwell's default google project the **google_project** workflow option must be set. In the similarly likely case that this service account can not access Cromwell's default google bucket, the **jes_gcs_root** workflow option should be set appropriately.


### Docker

It is possible to reference private docker images in DockerHub to be run on JES.
However, in order for the image to be pulled, the docker credentials with access to this image must be provided in the configuration file.


```
backend {
  default = "JES"
  providers {
    JES {
      actor-factory = "cromwell.backend.impl.local.LocalBackendLifecycleActorFactory"
      config {
        dockerhub {
          account = "mydockeraccount@mail.com"
          token = "mydockertoken"
        }
      }
    }
  }
}
```

It is now possible to reference an image only this account has access to:

```
task mytask {
  command {
    ...
  }
  runtime {
    docker: "private_repo/image"
    memory: "8 GB"
    cpu: "1"
  }
  ...
}
```

Note that if the docker image to be used is public there is no need to add this configuration.

### Monitoring

In order to monitor metrics (CPU, Memory, Disk usage...) about the VM during Call Runtime, a workflow option can be used to specify the path to a script that will run in the background and write its output to a log file.

```
{
  "monitoring_script": "gs://cromwell/monitoring/script.sh"
}
```

The output of this script will be written to a `monitoring.log` file that will be available in the call gcs bucket when the call completes.  This feature is meant to run a script in the background during long-running processes.  It's possible that if the task is very short that the log file does not flush before de-localization happens and you will end up with a zero byte file.

## GA4GH TES Backend
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

### Configuring
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

### Supported File Systems
Currently this backend only works with files on a Local or Shared File System. 

### Docker
This backend supports the following optional runtime attributes / workflow options for working with Docker:
* docker: Docker image to use such as "Ubuntu".
* dockerWorkingDir: defines the working directory in the container.

### CPU, Memory and Disk
This backend supports CPU, memory and disk size configuration through the use of the following runtime attributes / workflow options:
* cpu: defines the amount of CPU to use. Type: Integer. Ex: 4.
* memory: defines the amount of memory to use. Type: String. Ex: "4 GB" or "4096 MB"
* disk: defines the amount of disk to use. Type: String. Ex: "1 GB" or "1024 MB"

If they are not set, the TES backend may use default values.

## Sun GridEngine Backend

The GridEngine and similar backends use programs such as `qsub` to launch a job and will poll the filesystem to determine if a job is completed.

The backend is specified via the actor factory `ConfigBackendLifecycleActorFactory`:

```
backend {
  providers {
    SGE {
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

The job is launched with a configurable command such as:

```bash
qsub \
    -terse \
    -V \
    -b n \
    -N ${job_name} \
    -wd ${cwd} \
    -o ${out} \
    -e ${err} \
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
            -o ${out} \
            -e ${err} \
            ${script}
        """
      }
    }
  }
}
```

If the backend supports docker, another optional configuration key `backend.providers.<backend>.config.submit-docker` may be specified. When the WDL contains a docker runtime attribute, this command will be provided three additional variables:

* `docker` - The docker image name.
* `docker_cwd` - The path where `cwd` should be mounted within the docker container.
* `docker_cid` - The host path to which the [container ID file](https://docs.docker.com/engine/reference/run/#pid-equivalent) should be written.

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
            -o ${out} \
            -e ${err} \
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

## HtCondor Backend

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

```
chmod 755 ${script}
cat > ${cwd}/execution/dockerScript <<EOF
#!/bin/bash
docker run --rm -i -v ${cwd}:${docker_cwd} ${docker} /bin/bash ${script}
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

### Native Specifications
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

nativeSpecs attribute needs to be specified as String.

## Spark Backend

This backend adds support for execution of spark jobs in a workflow.

It supports the following Spark deploy modes:

*  Client deploy mode using the spark standalone cluster manager
*  Cluster deploy mode using the spark standalone cluster manager
*  Client deploy mode using Yarn resource manager
*  Cluster deploy mode using Yarn resource manager

### Configuring Spark Project

Cromwell's default configuration file is located at `core/src/main/resources/reference.conf`

To customize configuration it is recommended that one copies relevant stanzas from `core/src/main/resources/reference.conf` into a new file, modify it as appropriate, then pass it to Cromwell via:

java -Dconfig.file=/path/to/yourOverrides.conf cromwell.jar

Spark configuration stanza is as follows: 

```conf
Spark {
       actor-factory = "cromwell.backend.impl.spark.SparkBackendFactory"
       config {
         # Root directory where Cromwell writes job results.  This directory must be
         # visible and writeable by the Cromwell process as well as the jobs that Cromwell
         # launches.
         root: "cromwell-executions"

         filesystems {
           local {
             localization: [
               "hard-link", "soft-link", "copy"
             ]
           }
         }
		master: "local"
		deployMode: "client"
        }

      }
```
and add backend provider as Spark. 

```
backend {
  default = "Spark"
  providers {
  ....
```

### Configuring Spark Master and Deploy Mode

Default configuration is as follows:

```conf
Spark {
		......
		master: "local"
		deployMode: "client"

      }
```

However to use Spark in standalone cluster mode change `master: spark://hostname:6066` and `deployMode: cluster` similarly, for yarn change `master: yarn` and `deployMode: cluster` or `deployMode: client` to run in cluster or client mode respectively. 

### Spark runtime attributes

Supported runtime attributes for a Spark Job is as follows:
	
* executorCores (default value is 1)
* executorMemory (default value is "1GB", `Unit in MB or GB or TB.. ` )
* appMainClass ( Spark app/job entry point)
* numberOfExecutors ( Specific to cluster deploy mode)
* additionalArgs ( i.e to add additional configuration or parameters to spark-submit)

Sample usage:

```wdl
task sparkjob_with_yarn_cluster {
        .....
        
        runtime {
                appMainClass: "${entry_point}"
                executorMemory: "4GB"
                executorCores: "2"
                additionalArgs: "--conf '-Dsamjdk.compression_level=1 -Dsnappy.disable=true' ...."
        }
        
        .....
	}
```

### Spark Environment

The Spark backend assumes Spark is already installed, and it constructs the spark submit command with the `SPARK_HOME` environment variable if set. Otherwise backend creates command `spark-submit` without a fully qualified path to `spark-submit`.Also, it is important to set environment variable `HOSTNAME` to master machine ip or hostname, that is accessible by spark backend. That can be done by setting either in `~/.bashrc or profile like "export HOSTNAME=<machine ip>" `

Supported File Systems as follows: 

* Local File System
* Network File System
* Distributed file system

### Sample WDL
Next, create a WDL, and its json input like so:

```wdl
task sparkjob_with_yarn_cluster {
        File input_jar
        String input_1
        String output_base
        String entry_point
        Int cores
        String memory

        command {
                ${input_jar} ${input_1} ${output_base}
        }

        runtime {
                appMainClass: "${entry_point}"
                executorMemory: "${memory}"
                executorCores: "${cores}"
        }
        output {
                File out = "${output_base}"
          }
	}
```

and its accompanying json input as:

```json
{
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.memory": "4G",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.output_base":"/mnt/lustre/hadoop/home/yarn_cluster_output",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.entry_point": "com.org.spark.poc.nfs.SparkVowelLine",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.cores": "12",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.input_1": "/mnt/lustre/hadoop/home/inputfiles/sample.txt",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.input_jar": "/mnt/lustre/hadoop/home/inputjars/spark_hdfs.jar"
}
```