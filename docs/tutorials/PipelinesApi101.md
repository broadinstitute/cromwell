## Getting started on Google Cloud with the Genomics Pipelines API

## Pipelines API v2

### Basic Information

Initial support for Google [Pipelines API version 2](https://cloud.google.com/genomics/reference/rest/) was added in Cromwell 32.
Expect feature parity with v1 except:

* PAPI v2 private Docker support is equivalent to that in PAPI v1 but the configuration differs, please see [Docker configuration](http://cromwell.readthedocs.io/en/develop/filesystems/Google#Docker) for more details.
* The "refresh token" authentication mode is **NOT** supported on PAPI V2.

In addition, the following changes are to be expected:

* Error messages for failed jobs might differ from V1
* The Pipelines API log file content might differ from V1

### Setting up PAPIv2

For now the easiest way to try PAPIv2 is to migrate an existing set up from PAPIv1 (see below). After that, copy the PAPIv2 sample configuration in [cromwell.examples.conf](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.example.backends/PAPIv2.conf) in place of the PAPIv1 backend.

#### Permissions:

Google recommends using a service account to authenticate to GCP.  

You may create a service account using the `gcloud` command, consider running the following script and replace MY-GOOGLE-PROJECT:

```
#!/bin/bash
RANDOM_BUCKET_NAME=$(head /dev/urandom | tr -dc a-z | head -c 32 ; echo '')

#Create a new service account called "MyServiceAccount", and from the output of the command, take the email address that was generated
EMAIL=$(gcloud beta iam service-accounts create MyServiceAccount --description "to run cromwell"  --display-name "cromwell service account" --format json | jq '.email' | sed -e 's/\"//g')

# add all the roles to the service account
for i in storage.objectCreator storage.objectViewer genomics.pipelinesRunner genomics.admin iam.serviceAccountUser storage.objects.create
do
    gcloud projects add-iam-policy-binding MY-GOOGLE-PROJECT --member serviceAccount:"$EMAIL" --role roles/$i
done

# create a bucket to keep the execution directory
gsutil mb gs://"$RANDOM_BUCKET_NAME"

# give the service account write access to the new bucket
gsutil acl ch -u "$EMAIL":W gs://"$RANDOM_BUCKET_NAME"

# create a file that represents your service account.  KEEP THIS A SECRET.
gcloud iam service-accounts keys create sa.json --iam-account "$EMAIL"
```

## Pipelines API v1

### Deprecation warning

Please note that Google intends to deprecate PAPIv1 in the near future (circa mid 2019 or perhaps earlier). 

### Prerequisites

This tutorial page relies on completing the previous tutorial:

* [Downloading Prerequisites](FiveMinuteIntro.md)

### Goals

At the end of this tutorial you'll have run your first workflow against the Google Pipelines API.

### Let's get started!


**Configuring a Google Project**

Install the <a href="https://cloud.google.com/sdk/downloads" target="_blank">Google Cloud SDK</a>. 
Create a <a href="https://cloud.google.com/resource-manager/docs/creating-managing-projects" target="_blank">Google Cloud Project</a> and give it a project id (e.g. sample-project).  We’ll refer to this as `<google-project-id>` and your user login (e.g. username@gmail.com) as `<google-user-id>`.  

On your Google project, open up the <a href="https://console.developers.google.com/apis/library" target="_blank">API Manager</a> and enable the following APIs:
        
* Google Compute Engine
* Google Cloud Storage
* Genomics API

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
  default = "JES"
  providers {
    JES {
      actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
      config {
        // Google project
        project = "<google-project-id>"

        // Base bucket for workflow executions
        root = "gs://<google-bucket-name>/cromwell-execution"

        // Polling for completion backs-off gradually for slower-running jobs.
        // This is the maximum polling interval (in seconds):
        maximum-polling-interval = 600

        // Optional Dockerhub Credentials. Can be used to access private docker images.
        dockerhub {
          // account = ""
          // token = ""
        }

        genomics {
          // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
          // Pipelines and manipulate auth JSONs.
          auth = "application-default"
          // Endpoint for APIs, no reason to change this unless directed by Google.
          endpoint-url = "https://genomics.googleapis.com/"
          // This allows you to use an alternative service account to launch jobs, by default uses default service account
          compute-service-account = "default"

          // Pipelines v2 only: specify the number of times localization and delocalization operations should be attempted
          // There is no logic to determine if the error was transient or not, everything is retried upon failure
          // Defaults to 3
          localization-attempts = 3
        }

        filesystems {
          gcs {
            // A reference to a potentially different auth for manipulating files via engine functions.
            auth = "application-default"
            project = "<google-billing-project-id>"
          }
        }
      }
    }
  }
}
```

**Run Workflow**

`java -Dconfig.file=google.conf -jar cromwell-29.jar run hello.wdl -i hello.inputs`

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
