package cromwell.engine.finalcall

import cromwell.engine._
import cromwell.engine.backend.ExecutionHandle
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
trait FinalCall extends Scope {
  def workflow: WorkflowDescriptor
  def poll(implicit ec: ExecutionContext, executionHandle: ExecutionHandle): Future[ExecutionHandle]
  def execute(implicit ec: ExecutionContext): Future[ExecutionHandle]

  def prerequisiteCallNames: Set[LocallyQualifiedName] = throw new UnsupportedOperationException("prerequisiteCallNames not supported for FinalCalls")

  /**
    * This is handled as a special case by the WorkflowActor. We don't have to list the specially here.
    */
  def prerequisiteScopes: Set[Scope] = Set.empty

  val parent: Option[Scope] = Option(this.rootWorkflow)
}
