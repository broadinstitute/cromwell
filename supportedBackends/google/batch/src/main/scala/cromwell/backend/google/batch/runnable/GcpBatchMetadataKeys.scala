package cromwell.backend.google.batch.runnable

object GcpBatchMetadataKeys {
  val GoogleProject = "gcpBatch:googleProject"
  val ExecutionBucket = "gcpBatch:executionBucket"
  val MonitoringScript = "gcpBatch:monitoringScript"
  val MachineType = "gcpBatch:machineType"
  val Zone = "gcpBatch:zone"
  val InstanceName = "gcpBatch:instanceName"
  val MonitoringLog = "monitoringLog"
}
