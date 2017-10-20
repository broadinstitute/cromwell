A backend represents a way to run the user's command specified in the `task` section.  Cromwell allows for backends conforming to
the Cromwell backend specification to be plugged into the Cromwell engine.  Additionally, backends are included with the
Cromwell distribution:

* **Local / GridEngine / LSF / etc.** - Run jobs as subprocesses or via a dispatcher.  Supports launching in Docker containers. Use `bash`, `qsub`, `bsub`, etc. to run scripts.
* **Google Cloud** - Launch jobs on Google Compute Engine through the Google Genomics Pipelines API.
* **GA4GH TES** - Launch jobs on servers that support the GA4GH Task Execution Schema (TES).
* **HtCondor** - Allows to execute jobs using HTCondor.
* **Spark** - Adds support for execution of spark jobs.

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
* [Google Backend](#google-jes-backend) uses the [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem)

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

On the Google Pipelines backend the GCS (Google Cloud Storage) filesystem is used for the root of the workflow execution.
On the Local, SGE, and associated backends any GCS URI will be downloaded locally.  For the Google backend the `jes_gcs_root` [workflow option](#workflow-options) will take
precedence over the `root` specified at `backend.providers.JES.config.root` in the configuration file. Google Cloud Storage URIs are the only acceptable values for `File` inputs for
workflows using the Google backend.