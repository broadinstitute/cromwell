package cromwell.engine


import wdl4s.Workflow

import scala.concurrent.{ExecutionContext, Future}

object CopyWorkflowOutputs {
  val Name: LocallyQualifiedName = "$final_call$copy_workflow_outputs"
}

/**
  * Final call implementation that copies workflow outputs to a specified destination.
  */
case class CopyWorkflowOutputs(workflow: WorkflowDescriptor) extends FinalCall {
  override def unqualifiedName = CopyWorkflowOutputs.Name
  override def rootWorkflow: Workflow = workflow.namespace.workflow
  override def execute(implicit ec: ExecutionContext): Future[Unit] = workflow.copyWorkflowOutputs
}
