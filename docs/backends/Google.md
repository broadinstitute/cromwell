## Google Cloud Backend

Google Genomics Pipelines API is a Docker-as-a-service from Google. It was formerly called JES (Job Execution Service) so you will see references to JES in the configuration files and code.

### Configuring Google Project

You'll need the following things to get started:

* A Google Project (Manage/create projects [here](https://console.developers.google.com/project))
* A Google Cloud Storage bucket (View/create buckets in your project [here](https://console.cloud.google.com/storage/browser))

On your Google project, open up the [API Manager](https://console.developers.google.com/apis/library) and enable the following APIs:

* Google Compute Engine
* Google Cloud Storage
* Genomics API

If your project is `my-project` your bucket is `gs://my-bucket/`, then update your [Cromwell configuration file](/configuring) as follows:

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

If your project has API quotas other than the defaults set the `genomics-api-queries-per-100-seconds` value to be the lesser of the `Queries per 100 seconds per user` and `Queries per 100 seconds` quotas. This value will be used to help tune Cromwell's rate of interaction with Pipelines API.

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
the `genomics` and `filesystems.gcs` sections within a Google configuration block must reference an auth defined in this block.
The auth for the `genomics` section governs the interactions with Google itself, while `filesystems.gcs` governs the localization
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

A **refresh_token** field must be specified in the [workflow options](/workflowoptions) when submitting the job.  Omitting this field will cause the workflow to fail.

The refresh token is passed to Google along with the `client-id` and `client-secret` pair specified in the corresponding entry in `auths`.

#### User Service Account

A [JSON key file for the service account](####service-acocunt) must be passed in via the **user_service_account_json** field in the [workflow options](/workflowoptions) when submitting the job. Omitting this field will cause the workflow to fail. The JSON should be passed as a string and will need to have no newlines and all instances of *"* and *\n* escaped. 

In the likely event that this service account does not have access to Cromwell's default google project the **google_project** workflow option must be set. In the similarly likely case that this service account can not access Cromwell's default google bucket, the **jes_gcs_root** workflow option should be set appropriately.


### Docker

It is possible to reference private docker images in DockerHub to be run on Pipelines API.
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
