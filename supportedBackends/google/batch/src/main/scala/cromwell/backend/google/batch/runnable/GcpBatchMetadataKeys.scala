package cromwell.backend.google.batch.runnable

object GcpBatchMetadataKeys {
  val GoogleProject = "gcpBatch:googleProject"
  val ExecutionBucket = "gcpBatch:executionBucket"
  val MonitoringScript = "gcpBatch:monitoringScript"
  val MonitoringLog = "monitoringLog"
  // TODO unused metadata keys below, this is probably not right and there may not be tests to cover these
  val MachineType = "gcpBatch:machineType"
  val Zone = "gcpBatch:zone"
  val InstanceName = "gcpBatch:instanceName"
}
