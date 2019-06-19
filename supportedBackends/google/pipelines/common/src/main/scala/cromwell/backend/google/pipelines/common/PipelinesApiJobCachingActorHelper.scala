package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.io.{PipelinesApiAttachedDisk, PipelinesApiWorkingDisk}
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path
import cromwell.services.metadata.CallMetadataKeys

import scala.language.postfixOps

trait PipelinesApiJobCachingActorHelper extends StandardCachingActorHelper {
  this: PipelinesApiAsyncBackendJobExecutionActor with JobLogging =>

  lazy val initializationData: PipelinesApiBackendInitializationData = {
    backendInitializationDataAs[PipelinesApiBackendInitializationData]
  }

  lazy val pipelinesConfiguration: PipelinesApiConfiguration = initializationData.papiConfiguration

  lazy val pipelinesApiCallPaths: PipelinesApiJobPaths = jobPaths.asInstanceOf[PipelinesApiJobPaths]

  lazy val runtimeAttributes = PipelinesApiRuntimeAttributes(validatedRuntimeAttributes, pipelinesConfiguration.runtimeConfig)

  lazy val workingDisk: PipelinesApiAttachedDisk = runtimeAttributes.disks.find(_.name == PipelinesApiWorkingDisk.Name).get

  lazy val callRootPath: Path = pipelinesApiCallPaths.callExecutionRoot
  lazy val returnCodeFilename: String = pipelinesApiCallPaths.returnCodeFilename
  lazy val returnCodeGcsPath: Path = pipelinesApiCallPaths.returnCode
  lazy val jesLogPath: Path = pipelinesApiCallPaths.jesLogPath

  lazy val maxPreemption: Int = runtimeAttributes.preemptible
  def preemptible: Boolean

  lazy val jesAttributes: PipelinesApiConfigurationAttributes = pipelinesConfiguration.papiAttributes

  lazy val defaultLabels: Labels = {
    val workflow = jobDescriptor.workflowDescriptor
    val call = jobDescriptor.taskCall
    val subWorkflow = workflow.callable
    val subWorkflowLabels = if (!subWorkflow.equals(workflow.rootWorkflow))
      Labels("cromwell-sub-workflow-name" -> subWorkflow.name)
    else
      Labels.empty

    val alias = call.localName
    val aliasLabels = if (!alias.equals(call.callable.name))
      Labels("wdl-call-alias" -> alias)
    else
      Labels.empty

    Labels(
      "cromwell-workflow-id" -> s"cromwell-${workflow.rootWorkflowId}",
      "wdl-task-name" -> call.callable.name
    ) ++ subWorkflowLabels ++ aliasLabels
  }

  lazy val originalLabels: Labels = defaultLabels

  lazy val backendLabels: Seq[GoogleLabel] = GoogleLabels.safeLabels(originalLabels.asTuple :_*)

  lazy val originalLabelEvents: Map[String, String] = originalLabels.value map { l => s"${CallMetadataKeys.Labels}:${l.key}" -> l.value } toMap

  override protected def nonStandardMetadata: Map[String, Any] = {
    val googleProject = initializationData
      .workflowPaths
      .workflowDescriptor
      .workflowOptions
      .get(WorkflowOptionKeys.GoogleProject)
      .getOrElse(jesAttributes.project)

    Map(
      PipelinesApiMetadataKeys.GoogleProject -> googleProject,
      PipelinesApiMetadataKeys.ExecutionBucket -> initializationData.workflowPaths.executionRootString,
      PipelinesApiMetadataKeys.EndpointUrl -> jesAttributes.endpointUrl,
      "preemptible" -> preemptible
    ) ++ originalLabelEvents
  }
}
