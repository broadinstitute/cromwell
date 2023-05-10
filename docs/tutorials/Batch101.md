## Getting started on Google Cloud with Batch

## Batch

### Basic Information

Batch is a fully managed service that lets you schedule, queue, and execute batch processing workloads on Google Cloud resources. 
Batch provisions resources and manages capacity on your behalf, allowing your batch workloads to run at scale.

### Setting up Batch

#### Permissions:

### Prerequisites

This tutorial page relies on completing the previous tutorial:

* [Downloading Prerequisites](FiveMinuteIntro.md)

### Goals

At the end of this tutorial you'll have run your first workflow against the Google Batch API.

### Let's get started!


**Configuring a Google Project**

Install the <a href="https://cloud.google.com/sdk/downloads" target="_blank">Google Cloud SDK</a>.
Create a <a href="https://cloud.google.com/resource-manager/docs/creating-managing-projects" target="_blank">Google Cloud Project</a> and give it a project id (e.g. sample-project).  We’ll refer to this as `<google-project-id>` and your user login (e.g. username@gmail.com) as `<google-user-id>`.

On your Google project, open up the <a href="https://console.developers.google.com/apis/library" target="_blank">API Manager</a> and enable the following APIs:

* Google Compute Engine API
* Cloud Storage
* Google Cloud Batch API

Authenticate to Google Cloud Platform  
`gcloud auth login <google-user-id>`

Set your default account (will require to login again)  
`gcloud auth application-default login`

Set your default project  
`gcloud config set project <google-project-id>`

Create a Google Cloud Storage (GCS) bucket to hold Cromwell execution directories.
We will refer to this bucket as `google-bucket-name`, and the full identifier as `gs://google-bucket-name`.  
`gsutil mb gs://<google-bucket-name>`


**Workflow Source Files**

Copy over the sample `hello.wdl` and `hello.inputs` files to the same directory as the Cromwell jar.
This workflow takes a string value as specified in the inputs file and writes it to stdout.


***hello.wdl***
```
task hello {
  String addressee  
  command {
    echo "Hello ${addressee}! Welcome to Cromwell . . . on Google Cloud!"  
  }
  output {
    String message = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow wf_hello {
  call hello

  output {
     hello.message
  }
}
```

***hello.inputs***
```
{
  "wf_hello.hello.addressee": "World"
}
```

**Google Configuration File**

Copy over the sample `google.conf` file utilizing <a href="https://developers.google.com/identity/protocols/application-default-credentials" target="_blank">Application Default credentials</a> to the same directory that contains your sample WDL, inputs and Cromwell jar.
Replace `<google-project-id>` and `<google-bucket-name>`in the configuration file with the project id and bucket name. Replace `<google-billing-project-id>` with the project id that has to be billed for the request (more information for Requester Pays can be found at:
<a href="https://cloud.google.com/storage/docs/requester-pays" target="_blank">Requester Pays</a>)

***google.conf***
```
include required(classpath("application"))

google {

  application-name = "cromwell"

  auths = [
    {
      name = "application-default"
      scheme = "application_default"
    }
  ]
}

engine {
  filesystems {
    gcs {
      auth = "application-default"
      project = "<google-billing-project-id>"
    }
  }
}

backend {
  default = batch

  providers {
    batch {
      actor-factory = "cromwell.backend.google.pipelines.batch.GcpBatchBackendLifecycleActorFactory"
      config {
        # Google project
        project = "my-cromwell-workflows"

        # Base bucket for workflow executions
        root = "gs://my-cromwell-workflows-bucket"

        # Polling for completion backs-off gradually for slower-running jobs.
        # This is the maximum polling interval (in seconds):
        maximum-polling-interval = 600

        # Optional Dockerhub Credentials. Can be used to access private docker images.
        dockerhub {
          # account = ""
          # token = ""
        }

        # Optional configuration to use high security network (Virtual Private Cloud) for running jobs.
        # See https://cromwell.readthedocs.io/en/stable/backends/Google/ for more details.
        # virtual-private-cloud {
        #  network-label-key = "network-key"
        #  auth = "application-default"
        # }

        # Global pipeline timeout
        # Defaults to 7 days; max 30 days
        # batch-timeout = 7 days

        genomics {
          # A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
          # Batch Jobs and manipulate auth JSONs.
          auth = "application-default"


          // alternative service account to use on the launched compute instance
          // NOTE: If combined with service account authorization, both that service account and this service account
          // must be able to read and write to the 'root' GCS path
          compute-service-account = "default"

          # Location to submit jobs to Batch and store job metadata.
          location = "us-central1"

          # Specifies the minimum file size for `gsutil cp` to use parallel composite uploads during delocalization.
          # Parallel composite uploads can result in a significant improvement in delocalization speed for large files
          # but may introduce complexities in downloading such files from GCS, please see
          # https://cloud.google.com/storage/docs/gsutil/commands/cp#parallel-composite-uploads for more information.
          #
          # If set to 0 parallel composite uploads are turned off. The default Cromwell configuration turns off
          # parallel composite uploads, this sample configuration turns it on for files of 150M or larger.
          parallel-composite-upload-threshold="150M"
        }

        filesystems {
          gcs {
            # A reference to a potentially different auth for manipulating files via engine functions.
            auth = "application-default"
            # Google project which will be billed for the requests
            project = "google-billing-project"

            caching {
              # When a cache hit is found, the following duplication strategy will be followed to use the cached outputs
              # Possible values: "copy", "reference". Defaults to "copy"
              # "copy": Copy the output files
              # "reference": DO NOT copy the output files but point to the original output files instead.
              #              Will still make sure than all the original output files exist and are accessible before
              #              going forward with the cache hit.
              duplication-strategy = "copy"
            }
          }
        }
      }
    }
  }
}
```

**Run Workflow**

`java -Dconfig.file=google.conf -jar cromwell-67.jar run hello.wdl -i hello.inputs`

**Outputs**

The end of your workflow logs should report the workflow outputs.

```
[info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "wf_hello.hello.message": "Hello World! Welcome to Cromwell . . . on Google Cloud!"
  },
  "id": "08213b40-bcf5-470d-b8b7-1d1a9dccb10e"
}
```

Success!

### Next steps

You might find the following tutorials interesting to tackle next:

* [Persisting Data Between Restarts](PersistentServer)
* [Server Mode](ServerMode.md)
