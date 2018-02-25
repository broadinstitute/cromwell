package cromwell.backend.impl.jes

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.io.JobPaths
import cromwell.core.path.Path
import cromwell.services.metadata.CallMetadataKeys

object JesJobPaths {
  val JesLogPathKey = "jesLog"
  val JesMonitoringKey = "monitoring"
  val JesExecParamName = "exec"
}

final case class JesJobPaths(override val workflowPaths: JesWorkflowPaths, jobKey: BackendJobDescriptorKey) extends JobPaths {

  def jesLogBasename = {
    val index = jobKey.index.map(s => s"-$s").getOrElse("")
    s"${jobKey.node.localName}$index"
  }

  override val returnCodeFilename: String = s"$jesLogBasename-rc.txt"
  // These and `jesLogBasename` above are `def`s rather than `val`s because they are referenced polymorphically from
  // the initialization code of the extended `JobPaths` trait, but this class will not have initialized its `val`s
  // at the time that code runs.
  override def defaultStdoutFilename: String = s"$jesLogBasename-stdout.log"
  override def defaultStderrFilename: String = s"$jesLogBasename-stderr.log"
  override val scriptFilename: String = s"${JesJobPaths.JesExecParamName}.sh"

  val jesLogFilename: String = s"$jesLogBasename.log"
  lazy val jesLogPath: Path = callExecutionRoot.resolve(jesLogFilename)

  val jesMonitoringLogFilename: String = s"${JesJobPaths.JesMonitoringKey}.log"
  lazy val jesMonitoringLogPath: Path = callExecutionRoot.resolve(jesMonitoringLogFilename)

  val jesMonitoringScriptFilename: String = s"${JesJobPaths.JesMonitoringKey}.sh"

  override lazy val customMetadataPaths = Map(
    CallMetadataKeys.BackendLogsPrefix + ":log" -> jesLogPath
  ) ++ (
    workflowPaths.monitoringScriptPath map { p => Map(JesMetadataKeys.MonitoringScript -> p,
                                                      JesMetadataKeys.MonitoringLog -> jesMonitoringLogPath) } getOrElse Map.empty
  )

  override lazy val customDetritusPaths: Map[String, Path] = Map(
    JesJobPaths.JesLogPathKey -> jesLogPath
  )

  override lazy val customLogPaths: Map[String, Path] = Map(
    JesJobPaths.JesLogPathKey -> jesLogPath
  )
}
