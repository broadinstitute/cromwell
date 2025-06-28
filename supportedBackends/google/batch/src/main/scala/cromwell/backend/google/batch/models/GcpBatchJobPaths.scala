package cromwell.backend.google.batch.models

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.google.batch.runnable.GcpBatchMetadataKeys
import cromwell.backend.google.batch.util.{MemoryRetryRunnable, MemoryRetryStandard}
import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.services.metadata.CallMetadataKeys

object GcpBatchJobPaths {
  val BatchLogPathKey = "log"
  val BatchMonitoringKey = "monitoring"
  val BatchMonitoringImageKey = "monitoringImage"
  val BatchExecParamName = "exec"
  val GcsTransferLibraryName = "gcs_transfer.sh"
  val GcsLocalizationScriptName = "gcs_localization.sh"
  val GcsDelocalizationScriptName = "gcs_delocalization.sh"
  val DrsLocalizationManifestName = "drs_manifest"

}
case class GcpBatchJobPaths(override val workflowPaths: GcpBatchWorkflowPaths,
                            jobKey: BackendJobDescriptorKey,
                            override val isCallCacheCopyAttempt: Boolean = false
) extends JobPaths {

  def batchLogBasename = {
    val index = jobKey.index
      .map(s => s"-$s")
      .getOrElse("")
    s"${jobKey.node.localName}$index"
  }

  val batchLogFilename: String = s"$batchLogBasename.log"
  lazy val batchLogPath: Path = callExecutionRoot.resolve(batchLogFilename)

  val batchMonitoringLogFilename: String = s"${GcpBatchJobPaths.BatchMonitoringKey}.log"
  lazy val batchMonitoringLogPath: Path = callExecutionRoot.resolve(batchMonitoringLogFilename)

  val batchMonitoringScriptFilename: String = s"${GcpBatchJobPaths.BatchMonitoringKey}.sh"
  val batchMonitoringImageScriptFilename: String = s"${GcpBatchJobPaths.BatchMonitoringImageKey}.sh"

  val jesMonitoringLogFilename: String = s"${GcpBatchJobPaths.BatchMonitoringKey}.log"
  lazy val jesMonitoringLogPath: Path = callExecutionRoot.resolve(jesMonitoringLogFilename)

  // Return the path to the batch log file, if we are configured to create the log file.
  private lazy val maybeBatchLogPath: Option[Path] =
    workflowPaths.gcpBatchConfiguration.batchAttributes.logsPolicy match {
      case GcpBatchLogsPolicy.Path => Some(batchLogPath)
      case _ => None
    }

  override lazy val memoryRetryError: Option[Path] =
    workflowPaths.gcpBatchConfiguration.batchAttributes.memoryRetryCheckMode match {
      case MemoryRetryRunnable => None
      case MemoryRetryStandard => maybeBatchLogPath
    }

  override lazy val customMetadataPaths = {
    val backendLogsMetadata = maybeBatchLogPath map { p: Path =>
      Map(CallMetadataKeys.BackendLogsPrefix + ":log" -> p)
    } getOrElse Map.empty

    val monitoringScriptMetadata = workflowPaths.monitoringScriptPath map { p =>
      Map(GcpBatchMetadataKeys.MonitoringScript -> p, GcpBatchMetadataKeys.MonitoringLog -> jesMonitoringLogPath)
    } getOrElse Map.empty

    backendLogsMetadata ++ monitoringScriptMetadata
  }

  // This is required to include this log file in the collection of those copied with cache hits.
  // Only include it here if we're configured to
  override lazy val customDetritusPaths: Map[String, Path] = maybeBatchLogPath map { p: Path =>
    Map(GcpBatchJobPaths.BatchLogPathKey -> p)
  } getOrElse Map.empty

  override lazy val customLogPaths: Map[String, Path] = Map(
    GcpBatchJobPaths.BatchLogPathKey -> batchLogPath
  )

  override def standardOutputAndErrorPaths: Map[String, Path] =
    super.standardOutputAndErrorPaths map { case (k, v) =>
      val updated = workflowPaths.standardStreamNameToFileNameMetadataMapper(this, k)
      k -> v.parent.resolve(updated)
    }

  override def forCallCacheCopyAttempts: JobPaths = this.copy(isCallCacheCopyAttempt = true)
}
