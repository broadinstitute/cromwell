package cromwell.services.metadata

object CallMetadataKeys {
  val ExecutionEvents = "executionEvents"
  val ExecutionStatus = "executionStatus"
  val RuntimeAttributes = "runtimeAttributes"
  val Inputs = "inputs"
  val Outputs = "outputs"
  val ReturnCode = "returnCode"
  val Backend = "backend"
  val Start = "start"
  val End = "end"
  val RetryableFailure = "retryableFailure"
  val Failures = "failures"
  val Stdout = "stdout"
  val Stderr = "stderr"
  val BackendLogsPrefix = "backendLogs"
  val BackendStatus = "backendStatus"
  val JobId = "jobId"
  val CallRoot = "callRoot"
  val SubWorkflowId = "subWorkflowId"
  val SubWorkflowMetadata = "subWorkflowMetadata"
}
