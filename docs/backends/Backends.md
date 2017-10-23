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

* [Local Backend](/local) and associated backends primarily use the [Shared Local Filesystem](SharedFilesystem).
* [Google Backend](Google) uses the [Google Cloud Storage Filesystem](google).

Note that while Local, SGE, LSF, etc. backends use the local or network filesystem for the directory structure of a workflow, they are able to localize inputs
from GCS paths if configured to use a GCS filesystem.  See [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem) for more details.