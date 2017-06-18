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
  val CallCaching = "callCaching"
  val Labels = "labels"

  object CallCachingKeys {
    val EffectiveModeKey = CallCaching + MetadataKey.KeySeparator + "effectiveCallCachingMode"
    val ReadResultMetadataKey = CallCaching + MetadataKey.KeySeparator + "result"
    val HitResultMetadataKey = CallCaching + MetadataKey.KeySeparator + "hit"
    val AllowReuseMetadataKey = CallCaching + MetadataKey.KeySeparator + "allowResultReuse"
    val HitFailuresKey = CallCaching + MetadataKey.KeySeparator + "hitFailures"
    val HashFailuresKey = CallCaching + MetadataKey.KeySeparator + "hashFailures"
    val HashesKey = CallCaching + MetadataKey.KeySeparator + "hashes"
  }
}
