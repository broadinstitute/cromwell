package cromwell.backend.google.pipelines.common

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.services.metadata.CallMetadataKeys

object PipelinesApiJobPaths {
  val JesLogPathKey = "jesLog"
  val JesMonitoringKey = "monitoring"
  val JesExecParamName = "exec"
  val GcsTransferLibraryName = "gcs_transfer.sh"
  val GcsLocalizationScriptName = "gcs_localization.sh"
  val GcsDelocalizationScriptName = "gcs_delocalization.sh"
}

final case class PipelinesApiJobPaths(override val workflowPaths: PipelinesApiWorkflowPaths, jobKey: BackendJobDescriptorKey) extends JobPaths {

  // `jesLogBasename` is a `def` rather than a `val` because it is referenced polymorphically from
  // the initialization code of the extended `JobPaths` trait, but this class will not have initialized its `val`s
  // at the time that code runs.
  def jesLogBasename = {
    val index = jobKey.index.map(s => s"-$s").getOrElse("")
    s"${jobKey.node.localName}$index"
  }

  val jesLogFilename: String = s"$jesLogBasename.log"
  lazy val jesLogPath: Path = callExecutionRoot.resolve(jesLogFilename)

  val jesMonitoringLogFilename: String = s"${PipelinesApiJobPaths.JesMonitoringKey}.log"
  lazy val jesMonitoringLogPath: Path = callExecutionRoot.resolve(jesMonitoringLogFilename)

  val jesMonitoringScriptFilename: String = s"${PipelinesApiJobPaths.JesMonitoringKey}.sh"

  override lazy val customMetadataPaths = Map(
    CallMetadataKeys.BackendLogsPrefix + ":log" -> jesLogPath
  ) ++ (
    workflowPaths.monitoringScriptPath map { p => Map(PipelinesApiMetadataKeys.MonitoringScript -> p,
                                                      PipelinesApiMetadataKeys.MonitoringLog -> jesMonitoringLogPath) } getOrElse Map.empty
  )

  override lazy val customDetritusPaths: Map[String, Path] = Map(
    PipelinesApiJobPaths.JesLogPathKey -> jesLogPath
  )

  override lazy val customLogPaths: Map[String, Path] = Map(
    PipelinesApiJobPaths.JesLogPathKey -> jesLogPath
  )

  override def standardOutputAndErrorPaths: Map[String, Path] = {
    super.standardOutputAndErrorPaths map { case (k, v) =>
      val updated = workflowPaths.standardStreamNameToFileNameMetadataMapper(this, k)
      k -> v.parent.resolve(updated)
    }
  }
}
