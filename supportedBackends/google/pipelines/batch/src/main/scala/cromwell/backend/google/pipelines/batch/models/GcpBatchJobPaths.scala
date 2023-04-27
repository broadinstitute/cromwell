package cromwell.backend.google.pipelines.batch.models

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.io.JobPaths
import cromwell.core.path.Path

object GcpBatchJobPaths {

  val BatchLogPathKey = "jesLog"
  val BatchMonitoringKey = "monitoring"
  val BatchMonitoringImageKey = "monitoringImage"
  val BatchExecParamName = "exec"
  val GcsTransferLibraryName = "gcs_transfer.sh"
  val GcsLocalizationScriptName = "gcs_localization.sh"
  val GcsDelocalizationScriptName = "gcs_delocalization.sh"
  val DrsLocalizationManifestName = "drs_manifest"

}
case class GcpBatchJobPaths(override val workflowPaths: GcpBatchWorkflowPaths, jobKey: BackendJobDescriptorKey, override val isCallCacheCopyAttempt: Boolean = false) extends JobPaths {

  def batchLogBasename = {
    val index = jobKey
      .index
      .map(s => s"-$s")
      .getOrElse("")
    s"${
      jobKey
        .node
        .localName
    }$index"
  }

  val batchLogFilename: String = s"$batchLogBasename.log"
  lazy val batchLogPath: Path = callExecutionRoot.resolve(batchLogFilename)

  val batchMonitoringLogFilename: String = s"${GcpBatchJobPaths.BatchMonitoringKey}.log"
  lazy val batchMonitoringLogPath: Path = callExecutionRoot.resolve(batchMonitoringLogFilename)

  val batchMonitoringScriptFilename: String = s"${GcpBatchJobPaths.BatchMonitoringKey}.sh"
  val batchMonitoringImageScriptFilename: String = s"${GcpBatchJobPaths.BatchMonitoringImageKey}.sh"


  override val returnCodeFilename: String = s"$batchLogBasename-rc.txt"
  override def defaultStdoutFilename: String = s"$batchLogBasename-stdout.log"
  override def defaultStderrFilename: String = s"$batchLogBasename-stderr.log"

  override def forCallCacheCopyAttempts: JobPaths = this.copy(isCallCacheCopyAttempt = true)
}
