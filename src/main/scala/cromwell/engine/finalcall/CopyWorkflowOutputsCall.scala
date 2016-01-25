package cromwell.engine.finalcall

import cromwell.engine._
import cromwell.engine.backend.{SuccessfulFinalCallExecution, ExecutionResult, ExecutionHandle}
import cromwell.engine.workflow.{FinalCallKey, ExecutionStoreKey}
import wdl4s.Workflow

import scala.concurrent.{ExecutionContext, Future}

object CopyWorkflowOutputsCall {
  val UnqualifiedName: LocallyQualifiedName = "$final_call$copy_workflow_outputs"
  def name(workflowName: String): LocallyQualifiedName = workflowName + ".$final_call$copy_workflow_outputs"

  def isCopyWorkflowOutputsCallFqn(fqn: String) = fqn.endsWith(UnqualifiedName)
  def storeKey(workflow: WorkflowDescriptor): ExecutionStoreKey = FinalCallKey(CopyWorkflowOutputsCall(workflow))
}

/**
  * Final call implementation that copies workflow outputs to a specified destination.
  */
case class CopyWorkflowOutputsCall(override val workflow: WorkflowDescriptor) extends FinalCall {
  override def unqualifiedName = CopyWorkflowOutputsCall.name(workflow.name)
  override def rootWorkflow: Workflow = workflow.namespace.workflow
  override def execute(implicit ec: ExecutionContext): Future[ExecutionHandle] = workflow.copyWorkflowOutputs map {
    _ => CopyWorkflowOutputsHandle
  }
  override def poll(implicit ec: ExecutionContext, current: ExecutionHandle): Future[ExecutionHandle] = current match {
    case CopyWorkflowOutputsHandle => Future.successful(CopyWorkflowOutputsHandle)
    case _ => Future.failed(new IllegalStateException("CopyWorkflowOutputsCall can only be polled with a CopyWorkflowOutputsHandle."))
  }
}

case object CopyWorkflowOutputsHandle extends ExecutionHandle {
  override val isDone = true
  override val result = SuccessfulFinalCallExecution
}
