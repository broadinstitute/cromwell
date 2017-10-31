# Options
Keys | Possible Values | Description
--|--|--
<a name="jes_gcs_root">jes_gcs_root</a> | `string`  | Where outputs of the workflow will be written.  Expects this to be a GCS URL (e.g. `gs://my-bucket/workflows`).  If this is not set, this defaults to the value within `backend.jes.config.root` in the [Configuration](../Configuring).
<a name="google_compute_service_account">google_compute_service_account</a> | `string` | Alternate service account to use on the compute instance (e.g. `my-new-svcacct@my-google-project.iam.gserviceaccount.com`).  If this is not set, this defaults to the value within `backend.jes.config.genomics.compute-service-account` in the [Configuration](../Configuring) if specified or `default` otherwise.
<a name="google_project">google_project</a> | `string` |  Google project used to execute this workflow.
<a name="refresh_token">refresh_token</a> |`string` |   Only used if `localizeWithRefreshToken` is specified in the [Configuration](../Configuring).
<a name="auth_bucket">auth_bucket</a> |`string` |     A GCS URL that only Cromwell can write to.  The Cromwell account is determined by the `google.authScheme` (and the corresponding `google.userAuth` and `google.serviceAuth`). Defaults to the the value in [jes_gcs_root](#jes_gcs_root).
<a name="monitoring_script">monitoring_script</a> |`string` |   Specifies a GCS URL to a script that will be invoked prior to the user command being run.  For example, if the value for monitoring_script is `"gs://bucket/script.sh"`, it will be invoked as `./script.sh > monitoring.log &`.  The value `monitoring.log` file will be automatically de-localized.

# Example
```json
{
  "jes_gcs_root": "gs://my-bucket/workflows",
  "google_project": "my_google_project",
  "refresh_token": "1/Fjf8gfJr5fdfNf9dk26fdn23FDm4x",
  "google_compute_service_account": " my-new-svcacct@my-google-project.iam.gserviceaccount.com"
  "auth_bucket": "gs://my-auth-bucket/private",
  "monitoring_script": "gs://bucket/script.sh"
}
```


TODO

- link to exact options in Configuration file
