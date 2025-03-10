package cromwell.backend.google.batch.models

sealed abstract class GcpBatchExitCode(val code: Int) extends Product with Serializable

/**
 * Represents the possible exit codes from Batch.
 *
 * See: https://cloud.google.com/batch/docs/troubleshooting#reserved-exit-codes
 */
object GcpBatchExitCode {
  case object VMPreemption extends GcpBatchExitCode(50001)

  case object VMReportingTimeout extends GcpBatchExitCode(50002)

  case object VMRebootedDuringExecution extends GcpBatchExitCode(50003)

  case object VMAndTaskAreUnresponsive extends GcpBatchExitCode(50004)

  case object TaskRunsOverMaximumRuntime extends GcpBatchExitCode(50005)

  case object VMRecreatedDuringExecution extends GcpBatchExitCode(50006)

  val values: List[GcpBatchExitCode] = List(
    VMPreemption,
    VMReportingTimeout,
    VMRebootedDuringExecution,
    VMAndTaskAreUnresponsive,
    TaskRunsOverMaximumRuntime,
    VMRecreatedDuringExecution
  )

  // Special case non-error code - this does not correspond to any messages we get from Batch,
  // this is what we use for tasks that don't error out at the Batch layer.
  case object Success extends GcpBatchExitCode(0)

  def fromEventMessage(message: String): Option[GcpBatchExitCode] =
    values.find { target =>
      message.contains(s"exit code ${target.code}")
    }

}
