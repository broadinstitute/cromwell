package cromwell.backend.google.batch.models

sealed trait GcpBatchLogsPolicy extends Product with Serializable

object GcpBatchLogsPolicy {
  case object CloudLogging extends GcpBatchLogsPolicy
  case object Path extends GcpBatchLogsPolicy
}
