# Google Cloud Storage (GCS)

## Overview 

Cromwell supports workflows referencing objects stored in [Google Cloud Storage](https://cloud.google.com/storage/)
The Cromwell configuration for the GCS is as follow:

```hocon
filesystems {
  gcs {
    # A reference to a potentially different auth for manipulating files via engine functions.
    auth = "application-default"

    # Requester pays configuration. Remove this section to disable support of requester pays in Cromwell
    requester-pays {
      # Google project which will be billed for the requests
      project = "google-billing-project"

      # Time for which requester pays information about a bucket is being cached
      # Set to 0 to disable caching
      cache-ttl = "20 minutes"
    }

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
See the [Call Caching documentation](../CallCaching.md) for more information.

## Requester Pays

GCS has a feature called Requester Pays (RP). This section describes how Cromwell supports it and the consequences on cost. Please first read the [official documentation](https://cloud.google.com/storage/docs/requester-pays) if you're not already familiar with it.

Support for requester pays in Cromwell needs to be explicitly configured. It is **disabled** by default.
To enable RP in Cromwell, add the following configuration to any (and all) `filesystems.gcs` stanza you want to enable requester pays for:

```hocon
requester-pays {
  # Time for which requester pays information about a bucket is being cached
  # Set to 0 to disable caching
  cache-ttl = "20 minutes"
}
```

We will discuss this value in a second. First, it is important to say that when RP is enabled, Cromwell takes extra steps to:
 
1) Determine if a bucket has requester pays enabled.
Cromwell makes a request to GCS to know if a bucket has RP. Because this can be expansive when multiplied by the number of files a Cromwell server has to interact with, this information can be cached for a finite duration.
This duration can be configured using the `cache-ttl` field above. Although not recommended, caching can be disabled by setting the value to `0`. This will incur an additional request for every request to GCS.
     
2) If it does, make sure the appropriate project gets billed for accessing it
If a bucket has RP, requests to this bucket must be billed to a google project. Here is how Cromwell determines what project will be used:

- If a `google_project` was set in the [workflow options](../wf_options/Google.md) when the workflow was submitted, this value will be used
- Otherwise, the value of the `requester-pays.project` field in the `gcs` filesystem configuration will be used
- Otherwise, if the machine Cromwell runs on is authenticated to gcloud and a default project is set, this value will be used

**Important #1**: In order for a project to be billable to access a bucket with requester pays, the credentials used need to have the `serviceusage.services.use` permission on this project. 

**Important #2**: Pipelines API version 1 does **not** support buckets with requester pays, so while Cromwell itself might be able to access bucket with RP, jobs running on Pipelines API V1 with file inputs and / or outputs will **not** work.
For full requester pays support, use the [Pipelines API v2 Cromwell backend](https://github.com/broadinstitute/cromwell/blob/develop/CHANGELOG.md#pipelines-api-v2). 
