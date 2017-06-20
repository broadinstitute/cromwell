package cromwell.backend.impl.jes

import akka.actor.Actor
import cromwell.backend.impl.jes.io.{JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path
import cromwell.services.metadata.CallMetadataKeys

import scala.language.postfixOps

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
    jesCallPaths.workflowPaths.monitoringPath map { path =>
      JesFileInput(s"$MonitoringParamName-in", path.pathAsString,
        JesWorkingDisk.MountPoint.resolve(JesMonitoringScript), workingDisk)
    }
  }
  lazy val monitoringOutput: Option[JesFileOutput] = monitoringScript map { _ => JesFileOutput(s"$MonitoringParamName-out",
    defaultMonitoringOutputPath.pathAsString, JesMonitoringLogFile, workingDisk)
  }

  lazy val defaultLabels: Labels = {
    val workflow = jobDescriptor.workflowDescriptor
    val call = jobDescriptor.call
    val subWorkflow = workflow.workflow
    val subWorkflowLabels = if (!subWorkflow.equals(workflow.rootWorkflow))
      Labels("cromwell-sub-workflow-name" -> subWorkflow.unqualifiedName)
    else
      Labels.empty

    val alias = call.unqualifiedName
    val aliasLabels = if (!alias.equals(call.task.name))
      Labels("wdl-call-alias" -> alias)
    else
      Labels.empty

    Labels(
      "cromwell-workflow-id" -> s"cromwell-${workflow.rootWorkflowId}",
      "cromwell-workflow-name" -> workflow.rootWorkflow.unqualifiedName,
      "wdl-task-name" -> call.task.name
    ) ++ subWorkflowLabels ++ aliasLabels
  }

  lazy val backendLabels: Labels = defaultLabels ++ workflowDescriptor.customLabels

  lazy val backendLabelEvents: Map[String, String] = backendLabels.value map { l => s"${CallMetadataKeys.Labels}:${l.key}" -> l.value } toMap

  override protected def nonStandardMetadata: Map[String, Any] = {
    Map(
      JesMetadataKeys.GoogleProject -> jesAttributes.project,
      JesMetadataKeys.ExecutionBucket -> jesAttributes.executionBucket,
      JesMetadataKeys.EndpointUrl -> jesAttributes.endpointUrl,
      "preemptible" -> preemptible
    ) ++ backendLabelEvents
  }
}
