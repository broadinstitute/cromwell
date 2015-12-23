package cromwell.binding

import cromwell.engine.WorkflowDescriptor
import cromwell.engine.workflow.{ExecutionStoreKey, FinalCallKey}

import scala.concurrent.{ExecutionContext, Future}

object FinalCall {

  implicit class FinalCallString(val fqn: FullyQualifiedName) extends AnyVal {
    /** Does this FQN conform to a final call? */
    def isFinalCall = fqn startsWith "$final_call"
    /** Convert this FQN to a `FinalCall` or throw if it can't be converted. */
    def toStoreKey(workflow: WorkflowDescriptor): ExecutionStoreKey = {
      fqn match {
        case CopyWorkflowOutputs.Name => FinalCallKey(CopyWorkflowOutputs(workflow))
        case _ => throw new IllegalArgumentException(s"Unrecognized call key $fqn")
      }
    }
  }
}

/** Scope representing a "final call" that is inserted after all other workflow executions. */
trait FinalCall extends Scope {
  def execute(implicit ec: ExecutionContext): Future[Unit]

  def prerequisiteCallNames: Set[LocallyQualifiedName] = ???

  /**
    * Note: This implementation makes final calls dependent on all "real" scopes in the workflow.  For a system with
    * at most one final call this should be fine, but if there are > 1 final calls should they all be made runnable
    * at the same time (consuming multiple slots in the thread pool with potentially time consuming work) or run
    * serially?  This implementation starts all FinalCalls at the same time.
    */
  def prerequisiteScopes: Set[Scope] = rootWorkflow.children.toSet

  val parent: Option[Scope] = None
}
