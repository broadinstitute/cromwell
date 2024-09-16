package cromwell.backend.google.batch.models

sealed trait GcpBatchLogsPolicy extends Product with Serializable

object GcpBatchLogsPolicy {
  case object CloudLogging extends GcpBatchLogsPolicy
}
