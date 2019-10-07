
**Google Cloud Backend**

Google Genomics Pipelines API is a Docker-as-a-service from Google. It was formerly called JES (Job Execution Service);
you may see outdated references to the older JES terminology in Cromwell configuration files and code.

This section offers detailed configuration instructions for using Cromwell with the Pipelines API in all supported
authentication modes. Before reading futher in this section please see the
[Getting started on Google Pipelines API](../tutorials/PipelinesApi101) for instructions common to all authentication modes
and detailed instructions for the application default authentication scheme in particular.
The instructions below assume you have created a Google Cloud Storage bucket and a Google project enabled for the appropriate APIs.

**Configuring Authentication**

The `google` stanza in the Cromwell configuration file defines how to authenticate to Google.  There are four different
authentication schemes that might be used:

* `application_default` (default, recommended) - Use [application default](https://developers.google.com/identity/protocols/application-default-credentials) credentials.
* `service_account` - Use a specific service account and key file (in PEM format) to authenticate.
* `user_account` - Authenticate as a user.
* `refresh_token` - Authenticate each individual workflow using a refresh token supplied in the workflow options.
* `user_service_account` - Authenticate each individual workflow using service account credentials supplied in the workflow options.

The `auths` block in the `google` stanza defines the authentication schemes within a Cromwell deployment:

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

These authentication schemes can be referenced by name within other portions of the configuration file.  For example, both
the `genomics` and `filesystems.gcs` sections within a Google configuration block must reference an auth defined in this block.
The auth for the `genomics` section governs the interactions with Google itself, while `filesystems.gcs` governs the localization
of data into and out of GCE VMs.

**Application Default Credentials**

By default, application default credentials will be used.  Only `name` and `scheme` are required for application default credentials.

To authenticate, run the following commands from your command line (requires [gcloud](https://cloud.google.com/sdk/gcloud/)):

```
$ gcloud auth login
$ gcloud config set project my-project
```

**Service Account**

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

While technically not part of Service Account authentication mode, one can also override the default service account that the compute VM is started with via the configuration option `JES.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`.  The service account you provide must have been granted Service Account Actor role to Cromwell's primary service account. As this only affects Google Pipelines API and not GCS, it's important that this service account, and the service account specified in `JES.config.genomics.auth` can both read/write the location specified by `JES.config.root`

**Refresh Token**

A **refresh_token** field must be specified in the [Workflow Options](../wf_options/Google.md) when submitting the job.  Omitting this field will cause the workflow to fail.

The refresh token is passed to Google along with the `client-id` and `client-secret` pair specified in the corresponding entry in `auths`.

**User Service Account**

A [JSON key file for the service account](../wf_options/Google.md) must be passed in via the `user_service_account_json` field in the [Workflow Options](../wf_options/Google.md) when submitting the job. Omitting this field will cause the workflow to fail. The JSON should be passed as a string and will need to have no newlines and all instances of `"` and `\n` escaped. 

In the likely event that this service account does not have access to Cromwell's default google project the `google_project` workflow option must be set. In the similarly likely case that this service account can not access Cromwell's default google bucket, the `jes_gcs_root` workflow option should be set appropriately.

For information on the interaction of `user_service_account_json` with private Docker images please see the `Docker` section below.  

**Docker**

It's possible to reference private Docker images to which only particular Docker Hub accounts have access:

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

In order for a private image to be used the appropriate Docker configuration must be provided. If the Docker images being used
are public there is no need to add this configuration.

For Pipelines API (PAPI) version 1:
```
backend {
  default = "PAPIv1"
  providers {
    PAPIv1 {
      actor-factory = "cromwell.backend.google.pipelines.v1alpha2.PipelinesApiLifecycleActorFactory"
      config {
        dockerhub {
          token = "base64-encoded-docker-hub-username:password"
        }
      }
    }
  }
}
```

`token` is the standard base64-encoded username:password for the appropriate Docker Hub account.

For PAPI version 2:

```
backend {
  default = "PAPIv2"
  providers {
    PAPIv2 {
      actor-factory = "cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory"
      config {
        dockerhub {
          token = "base64-encoded-docker-hub-username:password"
          key-name = "name/of/the/kms/key/used/for/encrypting/and/decrypting/the/docker/hub/token"
          auth = "reference-to-the-auth-cromwell-should-use-for-kms-encryption"
        }
      }
    }
  }
}
```

`key-name` is the name of the Google KMS key Cromwell should use for encrypting the Docker `token` before including it
in the PAPI job execution request. This `key-name` will also be included in the PAPI job execution
request and will be used by PAPI to decrypt the Docker token used by `docker login` to enable access to the private Docker image.
 
`auth` is a reference to the name of an authorization in the `auths` block of Cromwell's `google` config.
Cromwell will use this authorization for encrypting the Google KMS key.

The equivalents of `key-name`, `token` and `auth` can also be specified in workflow options which take
precedence over values specified in configuration. The corresponding workflow options are named `docker_credentials_key_name`,
`docker_credentials_token`, and `user_service_account_json`. While the config value `auth` refers to an auth defined in the 
`google.auths` stanza elsewhere in Cromwell's
configuration, `user_service_account_json` is expected to be a literal escaped Google service account auth JSON.
See the `User Service Account` section above for more information on using user service accounts.
If the key, token or auth value is provided in workflow options then the corresponding private Docker configuration value
is not required, and vice versa. Also note that for the `user_service_account_json` workflow option to work an auth of type `user_service_account`
must be defined in Cromwell's `google.auths` stanza; more details in the `User Service Account` section above.

Example PAPI v2 workflow options for private Docker configuration:

```
{
  "docker_credentials_key_name": "name/of/the/kms/key/used/for/encrypting/and/decrypting/the/docker/hub/token",
  "docker_credentials_token": "base64_username:password",
  "user_service_account_json": "<properly escaped user service account JSON file>"
}
```

Important

If any of the three private Docker configuration values of key name, auth, or Docker token are missing, PAPI v2 will not perform a `docker login`.
If the Docker image to be pulled is not public the `docker pull` will fail which will cause the overall job to fail.

If using any of these private Docker workflow options it is advisable to add
them to the `workflow-options.encrypted-fields` list in Cromwell configuration.


**Monitoring**

In order to monitor metrics (CPU, Memory, Disk usage...) about the VM during Call Runtime, a workflow option can be used to specify the path to a script that will run in the background and write its output to a log file.

```
{
  "monitoring_script": "gs://cromwell/monitoring/script.sh"
}
```

The output of this script will be written to a `monitoring.log` file that will be available in the call gcs bucket when the call completes.  This feature is meant to run a script in the background during long-running processes.  It's possible that if the task is very short that the log file does not flush before de-localization happens and you will end up with a zero byte file.

**Google Cloud Storage Filesystem**

On the Google Pipelines backend the GCS (Google Cloud Storage) filesystem is used for the root of the workflow execution.
On the Local, SGE, and associated backends any GCS URI will be downloaded locally.  For the Google backend the `jes_gcs_root` [Workflow Option](../wf_options/Google) will take
precedence over the `root` specified at `backend.providers.JES.config.root` in the configuration file. Google Cloud Storage URIs are the only acceptable values for `File` inputs for
workflows using the Google backend.

**Retry with More Memory**

With `memory-retry` you can specify an array of strings which when encountered in the `stderr` file by Cromwell, allows the task to be retried with more memory.
The optional `multiplier` config specifies the factor by which the memory should be multiplied while retrying. This multiplier should be greater than 1.0. 
If the value is not mentioned in config, it will default to 2.0. The retry will be counted against the `maxRetries` count mentioned in the `runtimeAtrributes` in the task. 
For example,
```hocon
backend.providers.Papiv2.config {
  memory-retry {
    error-keys = ["OutOfMemoryError", "Killed"]
    multiplier = 1.1
  }
}
```  
this tells Cromwell to retry the task with 1.1x memory when it sees either `OutOfMemoryError` or `Killed` in the `stderr` file. If the task has 
runtime attributes as below 
```hocon
runtimeAtrributes {
  memory: "1 GB"
  continueOnReturnCode: true
  maxRetries: 1
}
``` 
the task will be retried at max 1 more time, and this time with "1.1 GB" memory. Please note that Pipelines API will adjust the memory value based on their
standards for memory for a VM. So it's possible that even though the request says 1.1 GB memory, it actually allocated a bit more memory to the VM.

Two environment variables called `${MEM_UNIT}` and `${MEM_SIZE}` are also available inside the command block of a task,
making it easy to retrieve the new value of memory on the machine.


#### Google Labels

Every call run on the Pipelines API backend is given certain labels by default, so that Google resources can be queried by these labels later. 
The current default label set automatically applied is:

| Key | Value | Example | Notes |
|-----|-------|---------|-------|
| cromwell-workflow-id | The Cromwell ID given to the root workflow (i.e. the ID returned by Cromwell on submission) | cromwell-d4b412c5-bf3d-4169-91b0-1b635ce47a26 | To fit the required [format](#label-format), we prefix with 'cromwell-' |
| cromwell-sub-workflow-name | The name of this job's sub-workflow | my-sub-workflow | Only present if the task is called in a subworkflow. |
| wdl-task-name | The name of the WDL task | my-task | |
| wdl-call-alias | The alias of the WDL call that created this job | my-task-1 | Only present if the task was called with an alias. |

Any custom labels provided as '`google_labels`' in the [workflow options](../wf_options/Google) are also applied to Google resources by the Pipelines API.

## Using NCBI Sequence Read Archive (SRA) Data

The v2alpha1 backend supports accessing [NCBI
SRA](https://www.ncbi.nlm.nih.gov/sra) accessions natively.  To configure this
support you'll need to enable it in your config file like so:

```hocon
filesystems {
  sra {
    class = "cromwell.filesystems.sra.SraPathBuilderFactory"
    docker-image = "fusera/fusera:alpine"
    ngc = "bmNiaV9nYXAfiwgAAAAAAAADBcHBDYQgEADAv1XQAGYXcfErUe5x0diCiFESA0Y8/VD8zTzrlXwMDEsoII9usPT5znZSmTqUohaSg5Gay14TbxsluMGOSBuqDEKefvbwCzv3BAAKoexb5uIbjjg7dq/p9mH7A5VTImxjAAAA"
  }
}
```

This filesystem has two required configuration options:
* `docker-image`: The [fusera](https://github.com/mitre/fusera) docker image to
  use to provide access.  This can be a custom image, but using the public
  [fusera/fusera:alpine](https://hub.docker.com/r/fusera/fusera/) image is
  recommended.
* `ngc`: A base-64 encoded NGC file.  This is provided through the NCBI
  interface.  Please see [the
  documentation](https://www.ncbi.nlm.nih.gov/books/NBK63512/#Download.are_downloaded_files_encrypted)
  for more information on obtaining your NGC.  The `ngc` value provided above
  is the sample credential file.

### Virtual Private Network

To run your jobs in a private network add the `virtual-private-cloud` stanza in the `config` stanza of the PAPI v2 backend:

```
backend {
  ...
  providers {
  	...
  	PapiV2 {
  	  actor-factory = "cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory"
  	  config {
  		...
  		virtual-private-cloud {
  	          network-label-key = "my-private-network"
  	          subnetwork-label-key = "my-private-subnetwork"
  	          auth = "reference-to-auth-scheme"
  	        }
  	    ...
  	  }  
      }
  }
}
```


The `network-label-key` and `subnetwork-label-key` should reference the keys in your project's labels whose value is the name of your private network
and subnetwork within that network respectively. `auth` should reference an auth scheme in the `google` stanza which will be used to get the project metadata from Google Cloud.
The `subnetwork-label-key` is an optional config.

For example, if your `virtual-private-cloud` config looks like the one above, and one of the labels in your project is

```
"my-private-network" = "vpc-network"
```

Cromwell will get labels from the project's metadata and look for a label whose key is `my-private-network`.
Then it will use the value of the label, which is `vpc-network` here, as the name of private network and run the jobs on this network.
If the network key is not present in the project's metadata Cromwell will fall back to running jobs on the default network.

### Parallel Composite Uploads

Cromwell can be configured to use GCS parallel composite uploads which can greatly improve delocalization performance. This feature
is turned off by default but can be enabled backend-wide by specifying a `gsutil`-compatible memory specification for the key
`genomics.parallel-composite-upload-threshold` in backend configuration. This memory value represents the minimum size an output file
must have to be a candidate for `gsutil` parallel composite uploading:

```
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory"
      config {
        ...
        genomics {
          ...
          parallel-composite-upload-threshold = 150M
          ...
        }
        ...
      }
    }
  }
}
```

Alternatively this threshold can be specified in workflow options using the key `parallel-composite-upload-threshold`,
which takes precedence over a setting in configuration. The default setting for this threshold is `0` which turns off
parallel composite uploads; a value of `0` can also be used in workflow options to turn off parallel composite uploads
in a Cromwell deployment where they are turned on in config.

#### Issues with composite files

Please see the [Google documentation](https://cloud.google.com/storage/docs/gsutil/commands/cp#parallel-composite-uploads)
describing the benefits and drawbacks of parallel composite uploads.

The actual error message observed when attempting to download a composite file on a system without a compiled `crcmod`
looks like the following:

```
/ # gsutil -o GSUtil:parallel_composite_upload_threshold=150M cp gs://my-bucket/composite.bam .
Copying gs://my-bucket/composite.bam...
==> NOTE: You are downloading one or more large file(s), which would
run significantly faster if you enabled sliced object downloads. This
feature is enabled by default but requires that compiled crcmod be
installed (see "gsutil help crcmod").

CommandException:
Downloading this composite object requires integrity checking with CRC32c,
but your crcmod installation isn't using the module's C extension, so the
hash computation will likely throttle download performance. For help
installing the extension, please see "gsutil help crcmod".

To download regardless of crcmod performance or to skip slow integrity
checks, see the "check_hashes" option in your boto config file.

NOTE: It is strongly recommended that you not disable integrity checks. Doing so
could allow data corruption to go undetected during uploading/downloading.
/ #
```

As the message states, the best option would be to have a compiled `crcmod` installed on the system.
Turning off integrity checks on downloads does get around this issue but really isn't a great idea.

#### Parallel composite uploads and call caching

Because the parallel composite upload threshold is not considered part of the hash used for call caching purposes, calls
which would be expected to generate non-composite outputs may call cache to results that did generate composite
outputs. Calls which are executed and not cached will always honor the parallel composite upload setting at the time of
their execution.
