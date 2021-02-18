**Google Cloud Backend**

Google Genomics Pipelines API is a Docker-as-a-service from Google. It was formerly called JES (Job Execution Service);
you may see outdated references to the older JES terminology in Cromwell configuration files and code.

This section offers detailed configuration instructions for using Cromwell with the Pipelines API in all supported
authentication modes. Before reading futher in this section please see the
[Getting started on Google Pipelines API](../tutorials/PipelinesApi101) for instructions common to all authentication modes
and detailed instructions for the application default authentication scheme in particular.
The instructions below assume you have created a Google Cloud Storage bucket and a Google project enabled for the appropriate APIs.

**Configuring Authentication**

The `google` stanza in the Cromwell configuration file defines how to authenticate to Google.  There are five different
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

***Application Default Credentials***

By default, application default credentials will be used.  Only `name` and `scheme` are required for application default credentials.

To authenticate, run the following commands from your command line (requires [gcloud](https://cloud.google.com/sdk/gcloud/)):

```
$ gcloud auth login
$ gcloud config set project my-project
```

***Service Account***

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

For PAPI version 2 alpha 1:

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

For PAPI version 2 beta:

```
backend {
  default = "PAPIv2"
  providers {
    PAPIv2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
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

**Pipeline timeout**

Google sets a default pipeline timeout of 7 days, after which the pipeline will abort. Setting `pipeline-timeout` overrides this limit to a maximum of 30 days.

```hocon
backend.providers.PAPIv2.config {
    pipeline-timeout: 14 days
}
```

**Enabling FUSE capabilities**

*This is a community contribution and not officially supported by the Cromwell team.*
By default Cromwell task containers doesn't allow to mount any FUSE filesystems. It happens because containers are launched without specific linux capabilities being enabled. 
Google pipelines backend supports running containers with the enabled capabilities and so does Cromwell. 

If you need to use fuses within task containers then you can set `enable_fuse` workflow option. 

```
{
    "enable_fuse": true
}
```

Differently you can enable support for fuses right in your backend configuration.

```
backend.providers.Papiv2.config {
    genomics {
        enable-fuse = true
    }
}
```

There is a list of limitations regarding the usage of FUSE filesystems:

+ Any inputs brought in via a FUSE filesystem will not be considered for call caching.
+ Any outputs stored via a FUSE filesystem will not be recreated if a task is replayed from a call-cache hit.
+ If the filesystem is writable, your job is potentially no longer idempotent - Cromwell may decide to retry your job for you, and you might get unforeseen file collisions or even incorrect results if that happens.

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

The v2alpha1 and v2beta backends support accessing [NCBI
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
  	  actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
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


### Custom Google Cloud SDK container
Cromwell can't use Google's container registry if VPC Perimeter is used in project.
Own repository can be used by adding `cloud-sdk-image-url` reference to used container:

```
google {
  ...
  cloud-sdk-image-url = "eu.gcr.io/your-project-id/cloudsdktool/cloud-sdk:275.0.0-slim"
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
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
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

### Migration from Google Cloud Genomics v2alpha1 to Google Cloud Life Sciences v2beta

1. If you currently run your workflows using Cloud Genomics v2alpha1 and would like to switch to Google Cloud Life 
Sciences v2beta, you will need to do a few changes to your configuration file: `actor-factory` value should be changed 
from `cromwell.backend.google.pipelines.v2alpha1.PipelinesApiLifecycleActorFactory` to `cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory`.
2. Parameter `genomics.endpoint-url` value should be changed from `https://genomics.googleapis.com/` to 
`https://lifesciences.googleapis.com/`.
3. Also you should add a new mandatory parameter `genomics.location` to your backend configuration. Currently Google Cloud 
Life Sciences API is available only in `us-central1` and `europe-west2` locations.

### Alpha support for WDL optional outputs on PAPI v2

Cromwell 53 adds alpha-quality support for WDL optional outputs on PAPI v2 backends. Constructs such as: 

```
  struct MyStruct {
    String name
    File? file
  }
  .
  .
  .
  output {
    File? file_does_not_exist = "does_not_exist"
    Pair[String, File?] pair_file_does_not_exist = ("this", "does_not_exist")
    Map[String, File?] map_file_does_not_exist = { "does_not_exist": "does_not_exist" }
    Array[File?] array_file_does_not_exist = ["does_not_exist"]
    MyStruct struct_file_does_not_exist = object { name: "this", file: "does_not_exist" } 
  }
```

will not produce errors if the file `does_not_exist` does not exist. This support for optional files is considered alpha
quality for two reasons:

1. As seen in the example above, support for optional files extends to complex WDL types but there is a restriction that
all `File` components of non-primitive types must be optional. e.g. Cromwell would not allow the assignment of a 
missing file to the right side of a pair of type `Pair[File, File?]` since the left member of the pair is a non-optional
file. This restriction exists solely due to technical limitations in how type evaluation works in Cromwell today and
may be removed in a future Cromwell release.

2. Call caching does not work for calls with empty optional outputs. Cromwell currently does not recognize
that it is okay for optional output files to be missing, will incorrectly claim that any cache hits with missing 
optional output files are unusable, and will proceed to search for more cache hits which if found will also be unusable,
before eventually giving up and running the job. This behavior may be corrected in a future Cromwell release.

### Reference Disk Support

Cromwell 55 and later support mounting reference disks from prebuilt GCP disk images as an alternative to localizing large
input reference files on PAPI v2. Please note the configuration of reference disk manifests has changed starting with
Cromwell 57 and now uses the format documented below. 

Within the `config` stanza of a PAPI v2 backend the `reference-disk-localization-manifests`
key specifies an array of reference disk manifests:  

```hocon
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
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

Using the first file in the manifest above as an example, assume a PAPI v2 backend is configured to use this manifest and the appropriate
`use_reference_disks` workflow option is set to `true` in the workflow submission. If a call in that workflow 
specifies the input `gs://my-references/enormous_reference.bam` and because that input matches the path of a file on the
reference image without the leading `gs://`, Cromwell would
arrange for a reference disk based on this image to be mounted and for the call's input to refer to the 
copy of the file on the reference disk, bypassing localization of the input.     

The Cromwell git repository includes a Java-based tool to facilitate the creation of manifests called
[CromwellRefdiskManifestCreatorApp](https://github.com/broadinstitute/cromwell/tree/develop/CromwellRefdiskManifestCreator).
Please see the help command of that tool for more details.

### Docker Image Cache Support

To optimize job execution time, Cromwell 55 and later support the use of Docker image caches on the PAPI v2 lifesciences beta backend. Docker image caches are not available on the PAPI v2 genomics alpha backend.
Configuration looks like:

```hocon
backend {
  ...
  providers {
    ...
    PapiV2 {
      actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
      config {
        ...
        docker-image-cache-manifest-file = "gs://path/to/a/docker/image/cache/manifest.json"
        ...
      }
    }
  }
}
```

Docker image cache manifest JSONs have a format like:

```json
{
  "biocontainers/samtools:1.3.1": "projects/broad-dsde-cromwell-dev/global/images/v1-docker-biocontainers-samtools-1-3-1",
  "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest": "projects/broad-dsde-cromwell-dev/global/images/v1-docker-gcr-io-gcp-runtimes-ubuntu-16-0-4-latest",
  ...
}
```

Docker image cache usage is an opt-in feature, so workflow submissions must specify this workflow option:

```json
{
  ...
  "use_docker_image_cache": true,
  ...
}
```

Individual tasks within a workflow can turn off Docker image caching through the use of a runtime attribute:

```wdl
task my_task {
  ...
  runtime {
    ...
    useDockerImageCache: false
  }
}
```

If Cromwell is running a workflow on PAPI v2 beta with Docker image caching enabled and a task specifies a
Docker image which corresponds to a configured Docker image cache JSON, Cromwell will arrange for the
job's VM to mount a disk built from the corresponding disk image. In the event that multiple
manifests describe disk images containing the specified Docker image, Cromwell will choose the disk image with the
smallest `diskSizeGb` value.

Conversely, Docker image caching can be turned off at the workflow level (either turned off explicitly or left at the
default setting of `false`) but turned on at the individual task level:

```wdl
task my_task {
  ...
  runtime {
    ...
    useDockerImageCache: true
  }
}
```

These settings could be useful for cost reasons: mounting Docker image caches adds nonzero cost
which might not be offset by eliminating Docker image pull times for long-running jobs.
