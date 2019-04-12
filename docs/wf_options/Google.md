# Google Pipelines API Workflow Options

These workflow options provide Google-specific information for workflows running tasks on the Google PAPI backend.

Keys | Possible Values | Description
--|--|--
`jes_gcs_root` | `string`  | Where outputs of the workflow will be written.  Expects this to be a GCS URL (e.g. `gs://my-bucket/workflows`).  If this is not set, this defaults to the value within `backend.jes.config.root` in the [Configuration](../Configuring).
`google_compute_service_account` | `string` | Alternate service account to use on the compute instance (e.g. `my-new-svcacct@my-google-project.iam.gserviceaccount.com`).  If this is not set, this defaults to the value within `backend.jes.config.genomics.compute-service-account` in the [Configuration](../Configuring) if specified or `default` otherwise.
`google_project` | `string` |  Google project used to execute this workflow.
`refresh_token` |`string` |   Only used if `localizeWithRefreshToken` is specified in the [Configuration](../Configuring).
`auth_bucket` |`string` |     A GCS URL that only Cromwell can write to.  The Cromwell account is determined by the `google.authScheme` (and the corresponding `google.userAuth` and `google.serviceAuth`). Defaults to the the value in [jes_gcs_root](#jes_gcs_root).
`monitoring_script` |`string` |   Specifies a GCS URL to a script that will be invoked prior to the user command being run.  For example, if the value for monitoring_script is `"gs://bucket/script.sh"`, it will be invoked as `./script.sh > monitoring.log &`.  The value `monitoring.log` file will be automatically de-localized.
`monitoring_image` |`string` |   Specifies a Docker image to monitor the task. This image will run concurrently with the task container, and provides an alternative mechanism to `monitoring_script` (the latter runs *inside* the task container). For example, one can use `quay.io/broadinstitute/cromwell-monitor`, which reports cpu/memory/disk utilization metrics to [Stackdriver](https://cloud.google.com/monitoring/).
`google_labels` | `object` | An object containing only string values. Represent custom labels to send with PAPI job requests. Per the PAPI specification, each key and value must conform to the regex `[a-z]([-a-z0-9]*[a-z0-9])?`.

# Example
```json
{
  "jes_gcs_root": "gs://my-bucket/workflows",
  "google_project": "my_google_project",
  "refresh_token": "1/Fjf8gfJr5fdfNf9dk26fdn23FDm4x",
  "google_compute_service_account": " my-new-svcacct@my-google-project.iam.gserviceaccount.com"
  "auth_bucket": "gs://my-auth-bucket/private",
  "monitoring_script": "gs://bucket/script.sh",
  "monitoring_image": "quay.io/broadinstitute/cromwell-monitor",
  "google_labels": {
    "custom-label": "custom-value"
  }
}
```
