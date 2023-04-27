package cromwell.backend.google.pipelines.batch.runnable

object GcpBatchApiMetadataKeys {
  val GoogleProject = "jes:googleProject"
  val ExecutionBucket = "jes:executionBucket"
  val EndpointUrl = "jes:endpointUrl"
  val MonitoringScript = "jes:monitoringScript"
  val MachineType = "jes:machineType"
  val Zone = "jes:zone"
  val InstanceName = "jes:instanceName"
  val MonitoringLog = "monitoringLog"
}
