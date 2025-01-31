package cromwell.backend.google.batch.actors

import cromwell.backend.google.batch.io.{GcpBatchAttachedDisk, GcpBatchWorkingDisk}
import cromwell.backend.google.batch.runnable.{GcpBatchMetadataKeys, WorkflowOptionKeys}
import cromwell.backend.google.batch.models._
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path
import cromwell.services.metadata.CallMetadataKeys

import scala.language.postfixOps

trait GcpBatchJobCachingActorHelper extends StandardCachingActorHelper {
  this: GcpBatchAsyncBackendJobExecutionActor with JobLogging =>

  lazy val initializationData: GcpBackendInitializationData =
    backendInitializationDataAs[GcpBackendInitializationData]
  lazy val batchConfiguration: GcpBatchConfiguration = initializationData.gcpBatchConfiguration

  lazy val gcpBatchCallPaths: GcpBatchJobPaths = jobPaths.asInstanceOf[GcpBatchJobPaths]

  lazy val runtimeAttributes = GcpBatchRuntimeAttributes(
    validatedRuntimeAttributes,
    batchConfiguration.runtimeConfig
  )

  lazy val maxPreemption: Int = runtimeAttributes.preemptible

  def preemptible: Boolean

  lazy val workingDisk: GcpBatchAttachedDisk = runtimeAttributes.disks.find(_.name == GcpBatchWorkingDisk.Name).get

  lazy val callRootPath: Path = gcpBatchCallPaths.callExecutionRoot
  lazy val returnCodeFilename: String = gcpBatchCallPaths.returnCodeFilename
  lazy val returnCodeGcsPath: Path = gcpBatchCallPaths.returnCode
  lazy val gcpBatchLogPath: Path = gcpBatchCallPaths.batchLogPath
  lazy val memoryRetryRCFilename: String = gcpBatchCallPaths.memoryRetryRCFilename
  lazy val memoryRetryRCGcsPath: Path = gcpBatchCallPaths.memoryRetryRC

  lazy val batchAttributes: GcpBatchConfigurationAttributes = batchConfiguration.batchAttributes

  lazy val defaultLabels: Labels = {
    val workflow = jobDescriptor.workflowDescriptor
    val call = jobDescriptor.taskCall
    val subWorkflow = workflow.callable
    val subWorkflowLabels =
      if (!subWorkflow.equals(workflow.rootWorkflow))
        Labels("cromwell-sub-workflow-name" -> subWorkflow.name)
      else
        Labels.empty

    val alias = call.localName
    val aliasLabels =
      if (!alias.equals(call.callable.name))
        Labels("wdl-call-alias" -> alias)
      else
        Labels.empty

    Labels(
      "cromwell-workflow-id" -> s"cromwell-${workflow.rootWorkflowId}",
      "wdl-task-name" -> call.callable.name
    ) ++ subWorkflowLabels ++ aliasLabels
  }

  lazy val originalLabels: Labels = defaultLabels

  lazy val backendLabels: Seq[GcpLabel] = GcpLabel.safeLabels(originalLabels.asTuple: _*)

  lazy val originalLabelEvents: Map[String, String] = originalLabels.value map { l =>
    s"${CallMetadataKeys.Labels}:${l.key}" -> l.value
  } toMap

  override protected def nonStandardMetadata: Map[String, Any] = {
    val googleProject = initializationData.workflowPaths.workflowDescriptor.workflowOptions
      .get(WorkflowOptionKeys.GoogleProject)
      .getOrElse(batchAttributes.project)

    Map[String, Any](
      GcpBatchMetadataKeys.GoogleProject -> googleProject,
      GcpBatchMetadataKeys.ExecutionBucket -> initializationData.workflowPaths.executionRootString,
      "preemptible" -> preemptible
    ) ++ originalLabelEvents
  }

}
