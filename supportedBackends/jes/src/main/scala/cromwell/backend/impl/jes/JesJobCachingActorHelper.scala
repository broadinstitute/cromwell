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

  lazy val maxPreemption: Int = runtimeAttributes.preemptible
  def preemptible: Boolean

  lazy val jesAttributes: JesAttributes = jesConfiguration.jesAttributes

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
      "wdl-task-name" -> call.task.name
    ) ++ subWorkflowLabels ++ aliasLabels
  }

  lazy val originalLabels: Labels = defaultLabels ++ workflowDescriptor.customLabels

  lazy val backendLabels: Labels = GoogleLabels.toLabels(originalLabels.asTuple :_*)

  lazy val originalLabelEvents: Map[String, String] = originalLabels.value map { l => s"${CallMetadataKeys.Labels}:${l.key}" -> l.value } toMap

  lazy val backendLabelEvents: Map[String, String] = backendLabels.value map { l => s"${CallMetadataKeys.BackendLabels}:${l.key}" -> l.value } toMap

  override protected def nonStandardMetadata: Map[String, Any] = {
    val googleProject = initializationData
      .workflowPaths
      .workflowDescriptor
      .workflowOptions
      .get(JesMetadataKeys.GoogleProject)
      .getOrElse(jesAttributes.project) 

    Map(
      JesMetadataKeys.GoogleProject -> googleProject,
      JesMetadataKeys.ExecutionBucket -> initializationData.workflowPaths.executionRootString,
      JesMetadataKeys.EndpointUrl -> jesAttributes.endpointUrl,
      "preemptible" -> preemptible
    ) ++ backendLabelEvents ++ originalLabelEvents
  }
}
