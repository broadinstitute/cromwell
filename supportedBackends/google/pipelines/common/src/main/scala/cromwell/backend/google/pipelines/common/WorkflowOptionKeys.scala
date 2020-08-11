package cromwell.backend.google.pipelines.common

object WorkflowOptionKeys {
  val MonitoringScript = "monitoring_script"
  val MonitoringImage = "monitoring_image"
  val MonitoringImageScript = "monitoring_image_script"
  val EnableSSHAccess = "enable_ssh_access"
  val GoogleProject = "google_project"
  val GoogleComputeServiceAccount = "google_compute_service_account"
  val EnableFuse = "enable_fuse"
  val GoogleLegacyMachineSelection = "google_legacy_machine_selection"
}
