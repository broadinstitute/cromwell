package cromwell.backend.google.pipelines.common

/**
  * Google cloud scopes that don't have constants defined elsewhere in Google Cloud Java API.
  */
object GoogleCloudScopes {
  /**
    * More restricted version of com.google.api.services.cloudkms.v1.CloudKMSScopes.CLOUD_PLATFORM
    * Could use that scope to keep things simple, but docs say to use a more restricted scope:
    *
    *   https://cloud.google.com/kms/docs/accessing-the-api#google_compute_engine
    *
    * For some reason this scope isn't listed as a constant under CloudKMSScopes.
    */
  val KmsScope = "https://www.googleapis.com/auth/cloudkms"

  /**
    * Scope to write metrics to Stackdriver Monitoring API.
    * Used by the monitoring action.
    *
    * For some reason we couldn't find this scope within Google libraries
    */
  val MonitoringWrite = "https://www.googleapis.com/auth/monitoring.write"
}
