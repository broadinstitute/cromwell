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
    * Note: This implementation makes final calls dependent on all "real" scopes in the workflow.  For a system with
    * at most one final call this should be fine, but if there are > 1 final calls should they all be made runnable
    * at the same time (consuming multiple slots in the thread pool with potentially time consuming work) or run
    * serially?  This implementation starts all FinalCalls at the same time.
    */
  def prerequisiteScopes: Set[Scope] = rootWorkflow.children.toSet

  val parent: Option[Scope] = Option(this.rootWorkflow)
}
