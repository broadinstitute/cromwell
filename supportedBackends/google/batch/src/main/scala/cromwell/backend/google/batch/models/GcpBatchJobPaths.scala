package cromwell.backend.google.batch.models

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.google.batch.runnable.GcpBatchMetadataKeys
import cromwell.backend.io.JobPaths
import cromwell.core.path.Path

object GcpBatchJobPaths {

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

  override lazy val customMetadataPaths =
    workflowPaths.monitoringScriptPath map { p =>
      Map(GcpBatchMetadataKeys.MonitoringScript -> p, GcpBatchMetadataKeys.MonitoringLog -> jesMonitoringLogPath)
    } getOrElse Map.empty

  override def standardOutputAndErrorPaths: Map[String, Path] =
    super.standardOutputAndErrorPaths map { case (k, v) =>
      val updated = workflowPaths.standardStreamNameToFileNameMetadataMapper(this, k)
      k -> v.parent.resolve(updated)
    }

  override def forCallCacheCopyAttempts: JobPaths = this.copy(isCallCacheCopyAttempt = true)
}
