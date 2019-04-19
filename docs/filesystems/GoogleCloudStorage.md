# Google Cloud Storage (GCS)

## Overview 

Cromwell supports workflows referencing objects stored in [Google Cloud Storage](https://cloud.google.com/storage/).
The Cromwell configuration for GCS is as follow:

```hocon
filesystems {
  gcs {
    # A reference to a potentially different auth for manipulating files via engine functions.
    auth = "application-default"

    # Google project which will be billed for requests on buckets with requester pays enabled
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
```

- The `auth` field refers to the authentication schema that should be used to authenticate requests. See [here](../backends/Google.md) for more info.
- The `project` field has to do with the Requester Pays feature (see below).
- The `caching.duplication-strategy` field determines how Cromwell should behave w.r.t output files when call is being cached. The default strategy `copy` is to copy the file to its new call location. As mentioned, `reference` will not copy the file and simply point the results to the existing location.
See the [Call Caching documentation](../cromwell_features/CallCaching.md) for more information.

## Requester Pays

GCS has a feature called Requester Pays (RP). This section describes how Cromwell supports it and the consequences on cost. Please first read the [official documentation](https://cloud.google.com/storage/docs/requester-pays) if you're not already familiar with it.

The billing project Cromwell uses to access a bucket with requester pays is determined as follows:

- If a `google_project` was set in the [workflow options](../wf_options/Google.md) when the workflow was submitted, this value is used
- Otherwise, the value of the `project` field in the `gcs` filesystem configuration is used
- Otherwise, if the machine Cromwell runs on is authenticated using gcloud and a default project is set, this value will be used

**Important Note #1**: In order for a project to be billable to access a bucket with requester pays, the credentials used need to have the `serviceusage.services.use` permission on this project. 

**Important Note #2**: Pipelines API version 1 does **not** support buckets with requester pays, so while Cromwell itself might be able to access bucket with RP, jobs running on Pipelines API V1 with file inputs and / or outputs will **not** work.
For full requester pays support, use the [Pipelines API v2 Cromwell backend](https://github.com/broadinstitute/cromwell/blob/develop/CHANGELOG.md#pipelines-api-v2). 

**Important Note #3**: Access to requester pays buckets from Cromwell is seamless, this also means that Cromwell will not report in the logs or metadata when it access a bucket with requester pays. It is the user's responsibility to be aware of the extra cost of running workflows access requester pays buckets.
