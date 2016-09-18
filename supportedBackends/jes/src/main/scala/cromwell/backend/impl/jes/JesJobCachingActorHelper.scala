package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{Actor, ActorRef}
import better.files._
import cromwell.backend.callcaching.JobCachingActorHelper
import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.WorkflowOptionKeys
import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.core.logging.JobLogging

import scala.language.postfixOps

trait JesJobCachingActorHelper extends JobCachingActorHelper {
  this: Actor with JobLogging =>

  val ExecParamName = "exec"
  val MonitoringParamName = "monitoring"

  val JesMonitoringScript = JesWorkingDisk.MountPoint.resolve("monitoring.sh")
  val JesMonitoringLogFile = JesWorkingDisk.MountPoint.resolve("monitoring.log")

  def jesConfiguration: JesConfiguration

  def initializationData: JesBackendInitializationData

  def serviceRegistryActor: ActorRef

  def getPath(str: String) = jesCallPaths.gcsFileSystem.getPath(str)

  override lazy val configurationDescriptor = jesConfiguration.configurationDescriptor

  lazy val jesCallPaths = initializationData.workflowPaths.toJesCallPaths(jobDescriptor.key)

  lazy val runtimeAttributes = JesRuntimeAttributes(jobDescriptor.runtimeAttributes, jobLogger)

  lazy val retryable = jobDescriptor.key.attempt <= runtimeAttributes.preemptible
  lazy val workingDisk: JesAttachedDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get

  lazy val callRootPath: Path = jesCallPaths.callRootPath
  lazy val returnCodeFilename = jesCallPaths.returnCodeFilename
  lazy val returnCodeGcsPath = jesCallPaths.returnCodePath
  lazy val jesStdoutFile = jesCallPaths.stdoutPath
  lazy val jesStderrFile = jesCallPaths.stderrPath
  lazy val jesLogFilename = jesCallPaths.jesLogFilename
  lazy val defaultMonitoringOutputPath = callRootPath.resolve(JesMonitoringLogFile)

  lazy val maxPreemption = runtimeAttributes.preemptible
  lazy val preemptible: Boolean = jobDescriptor.key.attempt <= maxPreemption

  lazy val jesAttributes = jesConfiguration.jesAttributes
  // TODO: Move monitoring paths to JesCallPaths
  lazy val monitoringScript: Option[JesInput] = {
    jobDescriptor.workflowDescriptor.workflowOptions.get(WorkflowOptionKeys.MonitoringScript) map { path =>
      JesFileInput(s"$MonitoringParamName-in", getPath(path).toString,
        JesWorkingDisk.MountPoint.resolve(JesMonitoringScript), workingDisk)
    } toOption
  }
  lazy val monitoringOutput = monitoringScript map { _ => JesFileOutput(s"$MonitoringParamName-out",
    defaultMonitoringOutputPath.toString, File(JesMonitoringLogFile).path, workingDisk)
  }

  lazy val metadataKeyValues: Map[String, Any] = {
    val runtimeAttributesMetadata: Map[String, Any] = runtimeAttributes.asMap map {
      case (key, value) => s"runtimeAttributes:$key" -> value
    }

    var fileMetadata: Map[String, Any] = jesCallPaths.metadataPaths
    if (monitoringOutput.nonEmpty) {
      // TODO: Move this to JesCallPaths
      fileMetadata += JesMetadataKeys.MonitoringLog -> monitoringOutput.get.gcs
    }

    val otherMetadata = Map(
      JesMetadataKeys.GoogleProject -> jesAttributes.project,
      JesMetadataKeys.ExecutionBucket -> jesAttributes.executionBucket,
      JesMetadataKeys.EndpointUrl -> jesAttributes.endpointUrl,
      "preemptible" -> preemptible,
      "cache:allowResultReuse" -> true
    )

    runtimeAttributesMetadata ++ fileMetadata ++ otherMetadata
  }
}
