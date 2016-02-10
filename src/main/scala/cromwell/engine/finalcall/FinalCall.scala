package cromwell.engine.finalcall

import cromwell.engine._
import cromwell.engine.workflow.ExecutionStoreKey
import wdl4s.Scope

import scala.concurrent.{ExecutionContext, Future}

object FinalCall {

  implicit class FinalCallString(val fqn: FullyQualifiedName) extends AnyVal {
    /** Does this FQN conform to a final call? */
    def isFinalCall = fqn contains "$final_call"

    def storeKey(workflow: WorkflowDescriptor): ExecutionStoreKey = {
      if (CopyWorkflowOutputsCall.isCopyWorkflowOutputsCallFqn(fqn)) { CopyWorkflowOutputsCall.storeKey(workflow) }
      else throw new IllegalArgumentException(s"Cannot convert $fqn into a FinalCall store key for $workflow")
    }
  }
}

/** Scope representing a "final call" that is inserted after all other workflow executions. */
  //TODO: Does it break if we extend from the Call. Copying files should be a function of the backends I think.
  // Either way..currently this flow is broken
trait FinalCall extends Scope {
  def workflow: WorkflowDescriptor
  def poll(implicit ec: ExecutionContext, executionHandle: Any): Future[Any]
  def execute(implicit ec: ExecutionContext): Future[Any]

  override def prerequisiteCallNames: Set[LocallyQualifiedName] = throw new UnsupportedOperationException("prerequisiteCallNames not supported for FinalCalls")

  /**
    * This is handled as a special case by the WorkflowActor. We don't have to list the specially here.
    */
  override def prerequisiteScopes: Set[Scope] = Set.empty

  override val parent: Option[Scope] = Option(this.rootWorkflow)
}
