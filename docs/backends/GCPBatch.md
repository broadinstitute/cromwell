**Google Cloud Backend**

[//]:
Google Cloud Batch is a fully managed service that lets you schedule, queue, and execute batch processing workloads on Google Cloud resources. Batch provisions resources and manages capacity on your behalf, allowing your batch workloads to run at scale.

This section offers detailed configuration instructions for using Cromwell with the Google Cloud Batch in all supported
authentication modes. Before reading further in this section please see the
[Getting started on Google Cloud Batch](../tutorials/Batch101) for instructions common to all authentication modes
and detailed instructions for the application default authentication scheme in particular.
The instructions below assume you have created a Google Cloud Storage bucket and a Google project enabled for the appropriate APIs.

**Configuring Authentication**

The `google` stanza in the Cromwell configuration file defines how to authenticate to Google.  There are four different
authentication schemes that might be used:

* `application_default` (default, recommended) - Use [application default](https://developers.google.com/identity/protocols/application-default-credentials) credentials.
* `service_account` - Use a specific service account and key file (in PEM format) to authenticate.
* `user_account` - Authenticate as a user.
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
the `GCPBATCH` and `filesystems.gcs` sections within a Google configuration block must reference an auth defined in this block.
The auth for the `GCPBATCH` section governs the interactions with Google itself, while `filesystems.gcs` governs the localization
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

While technically not part of Service Account authentication mode, one can also override the default service account that the compute VM is started with via the configuration option `GCPBATCH.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`.  The service account you provide must have been granted Service Account Actor role to Cromwell's primary service account. As this only affects Google Batch API and not GCS, it's important that this service account, and the service account specified in `GCPBATCH.config.genomics.auth` can both read/write the location specified by `GCPBATCH.config.root`

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

In order for a private image to be used, Docker Hub credentials must be provided. If the Docker images being used
are public there is no need to add this configuration.

For Batch

```
backend {
  default = GCPBATCH
  providers {
    GCPBATCH {
      actor-factory = "cromwell.backend.google.batch.GcpBatchBackendLifecycleActorFactory"
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

**Monitoring**

In order to monitor metrics (CPU, Memory, Disk usage...) about the VM during Call Runtime, a workflow option can be used to specify the path to a script that will run in the background and write its output to a log file.

```
{
  "monitoring_script": "gs://cromwell/monitoring/script.sh"
}
```

The output of this script will be written to a `monitoring.log` file that will be available in the call gcs bucket when the call completes.  This feature is meant to run a script in the background during long-running processes.  It's possible that if the task is very short that the log file does not flush before de-localization happens and you will end up with a zero byte file.

**Google Cloud Storage Filesystem**

On the Google Batch backend the GCS (Google Cloud Storage) filesystem is used for the root of the workflow execution.
On the Local, SGE, and associated backends any GCS URI will be downloaded locally.  For the Google backend the `jes_gcs_root` [Workflow Option](../wf_options/Google) will take
precedence over the `root` specified at `backend.providers.JES.config.root` in the configuration file. Google Cloud Storage URIs are the only acceptable values for `File` inputs for
workflows using the Google backend.

**Batch timeout**

Google sets a default pipeline timeout of 7 days, after which the pipeline will abort. Setting `batch-timeout` overrides this limit to a maximum of 30 days.

```hocon
backend.providers.GCPBATCH.config {
    batch-timeout: 14 days
}
```

#### Google Labels

Every call run on the GCP Batch backend is given certain labels by default, so that Google resources can be queried by these labels later. 
The current default label set automatically applied is:

| Key | Value | Example | Notes |
|-----|-------|---------|-------|
| cromwell-workflow-id | The Cromwell ID given to the root workflow (i.e. the ID returned by Cromwell on submission) | cromwell-d4b412c5-bf3d-4169-91b0-1b635ce47a26 | To fit the required [format](#label-format), we prefix with 'cromwell-' |
| cromwell-sub-workflow-name | The name of this job's sub-workflow | my-sub-workflow | Only present if the task is called in a subworkflow. |
| wdl-task-name | The name of the WDL task | my-task | |
| wdl-call-alias | The alias of the WDL call that created this job | my-task-1 | Only present if the task was called with an alias. |

Any custom labels provided as '`google_labels`' in the [workflow options](../wf_options/Google) are also applied to Google resources by GCP Batch.

### Virtual Private Network

Cromwell can arrange for jobs to run in specific GCP private networks via the `config.virtual-private-cloud` stanza of a Batch backend.
There are two ways of specifying private networks:

* [Literal network and subnetwork values](#virtual-private-network-via-literals) that will apply to all projects
* [Google project labels](#virtual-private-network-via-labels) whose values in a particular Google project will specify the network and subnetwork

#### Virtual Private Network via Literals

```hocon
backend {
  ...
  providers {
    ...
    GCPBATCH {
      actor-factory = "cromwell.backend.google.batch.GcpBatchLifecycleActorFactory"
      config {
        ...
        virtual-private-cloud {
          network-name = "vpc-network"
          subnetwork-name = "vpc-subnetwork"
        }
        ...
      }
    }
  }
}
```

The `network-name` and `subnetwork-name` should reference the name of your private network and subnetwork within that
network respectively. The `subnetwork-name` is an optional config.

For example, if your `virtual-private-cloud` config looks like the one above, then Cromwell will use the value of the
configuration key, which is `vpc-network` here, as the name of private network and run the jobs on this network.
If the network name is not present in the config Cromwell will fall back to trying to run jobs on the default network.

If the `network-name` or `subnetwork-name` values contain the string `${projectId}` then that value will be replaced
by Cromwell with the name of the project running GCP Batch.

If the `network-name` does not contain a `/` then it will be prefixed with `projects/${projectId}/global/networks/`.

Cromwell will then pass the network and subnetwork values to GCP Batch. See the documentation for
[GCP Batch](https://cloud.google.com/batch/docs/networking-overview)
for more information on the various formats accepted for `network` and `subnetwork`.

#### Virtual Private Network via Labels

```hocon
backend {
  ...
  providers {
    ...
    GCPBATCH {
      actor-factory = "cromwell.backend.google.batch.GcpBatchLifecycleActorFactory"
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
If the network key is not present in the project's metadata Cromwell will fall back to trying to run jobs using literal
network labels, and then fall back to running on the default network.

### Custom Google Cloud SDK container

Cromwell can't use Google's container registry if VPC Perimeter is used in project.
Own repository can be used by adding `cloud-sdk-image-url` reference to used container:

```
google {
  ...
  cloud-sdk-image-url = "eu.gcr.io/your-project-id/cloudsdktool/cloud-sdk:354.0.0-alpine"
  cloud-sdk-image-size-gb = 1
}
```

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
    GCPBATCH {
      actor-factory = "cromwell.backend.google.batch.GcpBatchLifecycleActorFactory"
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


### Migration from Google Cloud Life Sciences v2beta to Google Cloud Batch

1. If you currently run your workflows using Cloud Genomics v2beta and would like to switch to Google Cloud Batch, you will need to do a few changes to your configuration file: `actor-factory` value should be changed 
from `cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory` to `cromwell.backend.google.batch.GcpBatchLifecycleActorFactory`.

2. You will need to remove the parameter `genomics.endpoint-url` and generate a new config file.

3. Google Cloud Batch is now available in a variety of regions. Please see the [Batch Locations](https://cloud.google.com/batch/docs/locations) for a list of supported regions


### Reference Disk Support

Cromwell 55 and later support mounting reference disks from prebuilt GCP disk images as an alternative to localizing large
input reference files on Batch. Please note the configuration of reference disk manifests has changed starting with
Cromwell 57 and now uses the format documented below. 

Within the `config` stanza of a Batch backend the `reference-disk-localization-manifests`
key specifies an array of reference disk manifests:  

```hocon
backend {
  ...
  providers {
    ...
    GCPBATCH {
      actor-factory = "cromwell.backend.google.batch.GcpBatchLifecycleActorFactory"
      config {
        ...
        reference-disk-localization-manifests = [
          {
            "imageIdentifier" : "projects/broad-dsde-cromwell-dev/global/images/broad-references-disk-image",
            "diskSizeGb" : 500,
            "files" : [ {
              "path" : "gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.nhr",
              "crc32c" : 407769621
            }, {
              "path" : "gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.sa",
              "crc32c" : 1902048083
            },
            ...
          },
          ...
        ]
        ...
      }
    }
  }
}
```

Reference disk usage is an opt-in feature, so workflow submissions must specify this workflow option:

```json
{
  ...
  "use_reference_disks": true,
  ...
}
```

Using the first file in the manifest above as an example, assume a Batch backend is configured to use this manifest and the appropriate
`use_reference_disks` workflow option is set to `true` in the workflow submission. If a call in that workflow 
specifies the input `gs://my-references/enormous_reference.bam` and because that input matches the path of a file on the
reference image without the leading `gs://`, Cromwell would
arrange for a reference disk based on this image to be mounted and for the call's input to refer to the 
copy of the file on the reference disk, bypassing localization of the input.     

The Cromwell git repository includes a Java-based tool to facilitate the creation of manifests called
[CromwellRefdiskManifestCreatorApp](https://github.com/broadinstitute/cromwell/tree/develop/CromwellRefdiskManifestCreator).
Please see the help command of that tool for more details.

Alternatively for public data stored under `gs://gcp-public-data--broad-references` there exists a shell script to
extract reference data to a new disk and then convert that disk to a public image. For more information see
[create_images.sh](https://github.com/broadinstitute/cromwell/tree/develop/scripts/reference_disks/create_images.sh).

