package cromwell.backend.impl.jes

import akka.actor.Actor
import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path

trait JesJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>

  val ExecParamName = "exec"
  val MonitoringParamName = "monitoring"

  val JesMonitoringScript: Path = JesWorkingDisk.MountPoint.resolve("monitoring.sh")
  val JesMonitoringLogFile: Path = JesWorkingDisk.MountPoint.resolve("monitoring.log")

  lazy val initializationData: JesBackendInitializationData = {
    backendInitializationDataAs[JesBackendInitializationData]
  }

  lazy val jesConfiguration: JesConfiguration = initializationData.jesConfiguration

  lazy val jesCallPaths: JesJobPaths = jobPaths.asInstanceOf[JesJobPaths]

  lazy val runtimeAttributes = JesRuntimeAttributes(validatedRuntimeAttributes, jesConfiguration.runtimeConfig)

  lazy val workingDisk: JesAttachedDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get

  lazy val callRootPath: Path = jesCallPaths.callExecutionRoot
  lazy val returnCodeFilename: String = jesCallPaths.returnCodeFilename
  lazy val returnCodeGcsPath: Path = jesCallPaths.returnCode
  lazy val jesStdoutFile: Path = jesCallPaths.stdout
  lazy val jesStderrFile: Path = jesCallPaths.stderr
  lazy val jesLogFilename: String = jesCallPaths.jesLogFilename
  lazy val defaultMonitoringOutputPath: Path = callRootPath.resolve(JesMonitoringLogFile)

  lazy val maxPreemption: Int = runtimeAttributes.preemptible
  def preemptible: Boolean

  lazy val jesAttributes: JesAttributes = jesConfiguration.jesAttributes
  lazy val monitoringScript: Option[JesInput] = {
    jesCallPaths.monitoringPath map { path =>
      JesFileInput(s"$MonitoringParamName-in", path.pathAsString,
        JesWorkingDisk.MountPoint.resolve(JesMonitoringScript), workingDisk)
    }
  }
  lazy val monitoringOutput: Option[JesFileOutput] = monitoringScript map { _ => JesFileOutput(s"$MonitoringParamName-out",
    defaultMonitoringOutputPath.pathAsString, JesMonitoringLogFile, workingDisk)
  }

  override protected def nonStandardMetadata: Map[String, Any] = {
    Map(
      JesMetadataKeys.GoogleProject -> jesAttributes.project,
      JesMetadataKeys.ExecutionBucket -> jesAttributes.executionBucket,
      JesMetadataKeys.EndpointUrl -> jesAttributes.endpointUrl,
      "preemptible" -> preemptible
    )
  }
}
