package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{Actor, ActorRef}
import better.files._
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.callcaching.JobCachingActorHelper
import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.core.logging.JobLogging

import scala.util.Try

trait JesJobCachingActorHelper extends JobCachingActorHelper {
  this: Actor with JobLogging =>

  val ExecParamName = "exec"
  val MonitoringParamName = "monitoring"

  val JesMonitoringScript = JesWorkingDisk.MountPoint.resolve("monitoring.sh")
  val JesMonitoringLogFile = JesWorkingDisk.MountPoint.resolve("monitoring.log")

  def jesConfiguration: JesConfiguration

  def initializationData: JesBackendInitializationData

  def serviceRegistryActor: ActorRef
  
  def workflowDescriptor: BackendWorkflowDescriptor

  def getPath(str: String): Try[Path] = jesCallPaths.getPath(str)

  override lazy val configurationDescriptor = jesConfiguration.configurationDescriptor

  lazy val jesCallPaths = {
    val workflowPaths = if (workflowDescriptor.breadCrumbs.isEmpty) {
      initializationData.workflowPaths
    } else {
      new JesWorkflowPaths(workflowDescriptor, jesConfiguration)(context.system)
    }
    
    workflowPaths.toJobPaths(jobDescriptor.key)
  }

  lazy val runtimeAttributes = JesRuntimeAttributes(jobDescriptor.runtimeAttributes, jobLogger)

  lazy val retryable = jobDescriptor.key.attempt <= runtimeAttributes.preemptible
  lazy val workingDisk: JesAttachedDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get

  lazy val callRootPath: Path = jesCallPaths.callExecutionRoot
  lazy val returnCodeFilename = jesCallPaths.returnCodeFilename
  lazy val returnCodeGcsPath = jesCallPaths.returnCode
  lazy val jesStdoutFile = jesCallPaths.stdout
  lazy val jesStderrFile = jesCallPaths.stderr
  lazy val jesLogFilename = jesCallPaths.jesLogFilename
  lazy val defaultMonitoringOutputPath = callRootPath.resolve(JesMonitoringLogFile)

  lazy val maxPreemption = runtimeAttributes.preemptible
  lazy val preemptible: Boolean = jobDescriptor.key.attempt <= maxPreemption

  lazy val jesAttributes = jesConfiguration.jesAttributes
  lazy val monitoringScript: Option[JesInput] = {
    jesCallPaths.monitoringPath map { path =>
      JesFileInput(s"$MonitoringParamName-in", path.toUri.toString,
        JesWorkingDisk.MountPoint.resolve(JesMonitoringScript), workingDisk)
    }
  }
  lazy val monitoringOutput = monitoringScript map { _ => JesFileOutput(s"$MonitoringParamName-out",
    defaultMonitoringOutputPath.toString, File(JesMonitoringLogFile).path, workingDisk)
  }

  // Implements CacheHitDuplicating.metadataKeyValues
  lazy val metadataKeyValues: Map[String, Any] = {
    val runtimeAttributesMetadata: Map[String, Any] = runtimeAttributes.asMap map {
      case (key, value) => s"runtimeAttributes:$key" -> value
    }
    
    val otherMetadata: Map[String, Any] = Map(
      JesMetadataKeys.GoogleProject -> jesAttributes.project,
      JesMetadataKeys.ExecutionBucket -> jesAttributes.executionBucket,
      JesMetadataKeys.EndpointUrl -> jesAttributes.endpointUrl,
      "preemptible" -> preemptible,
      "cache:allowResultReuse" -> true
    )

    runtimeAttributesMetadata ++ jesCallPaths.metadataPaths ++ otherMetadata
  }
}
