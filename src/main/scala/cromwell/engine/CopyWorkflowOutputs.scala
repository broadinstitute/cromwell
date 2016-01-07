package cromwell.engine

import cromwell.binding._

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

  /** TODO replace snarky println with call to Thibault's code that actually copies files. */
  override def execute(implicit ec: ExecutionContext): Future[Unit] = Future.successful( { println("I have copied your outputs.  You're welcome.") } )
}
