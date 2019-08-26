package cromwell.services.metadata

object CallMetadataKeys {

  val CompressedDockerSize = "compressedDockerSize"
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
  val CallCaching = "callCaching"
  val BackendLabels = "backendLabels"
  val Labels = "labels"
  val CommandLine = "commandLine"

  object CallCachingKeys {
    val EffectiveModeKey = CallCaching + MetadataKey.KeySeparator + "effectiveCallCachingMode"
    val ReadResultMetadataKey = CallCaching + MetadataKey.KeySeparator + "result"
    val HitResultMetadataKey = CallCaching + MetadataKey.KeySeparator + "hit"
    val AllowReuseMetadataKey = CallCaching + MetadataKey.KeySeparator + "allowResultReuse"
    val HashFailuresKey = CallCaching + MetadataKey.KeySeparator + "hashFailures"
    val HashesKey = CallCaching + MetadataKey.KeySeparator + "hashes"
  }
}
