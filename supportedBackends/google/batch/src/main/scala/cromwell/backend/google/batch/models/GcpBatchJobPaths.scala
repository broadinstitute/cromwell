package cromwell.backend.google.batch.models

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.google.batch.runnable.GcpBatchMetadataKeys
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

  override lazy val customMetadataPaths = Map(
    CallMetadataKeys.BackendLogsPrefix + ":log" -> batchLogPath
  ) ++ (
    workflowPaths.monitoringScriptPath map { p =>
      Map(GcpBatchMetadataKeys.MonitoringScript -> p, GcpBatchMetadataKeys.MonitoringLog -> jesMonitoringLogPath)
    } getOrElse Map.empty
  )

  // This is required to include this log file in the collection of those copied with cache hits.
  // Only include it here if we're configured to
  override lazy val customDetritusPaths: Map[String, Path] =
    workflowPaths.gcpBatchConfiguration.batchAttributes.logsPolicy match {
      case GcpBatchLogsPolicy.Path =>
        Map(GcpBatchJobPaths.BatchLogPathKey -> batchLogPath)
      case _ => Map.empty
    }

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
