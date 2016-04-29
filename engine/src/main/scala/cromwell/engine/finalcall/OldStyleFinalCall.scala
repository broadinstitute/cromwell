package cromwell.engine.finalcall

import cromwell.core.LocallyQualifiedName
import cromwell.engine._
import cromwell.engine.backend.{OldStyleWorkflowDescriptor, OldStyleExecutionHandle, OldStyleSuccessfulFinalCallExecution}
import cromwell.engine.workflow.{ExecutionStoreKey, FinalCallKey}
import cromwell.webservice.WorkflowMetadataResponse
import wdl4s.{Scope, Workflow}

import scala.collection.immutable.Traversable
import scala.concurrent.{ExecutionContext, Future}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleFinalCall {

  private val finalCallCompanions: Traversable[FinalCallCompanion[_ <: OldStyleFinalCall]] = List(
    CopyWorkflowLogCall,
    CopyWorkflowOutputsCall,
    CopyCallLogsCall
  )

  def createFinalCalls(descriptor: OldStyleWorkflowDescriptor): Traversable[OldStyleFinalCall] = {
    finalCallCompanions.map(_.createCall(descriptor))
  }

  implicit class FinalCallString(val fqn: FullyQualifiedName) extends AnyVal {
    /** Does this FQN conform to a final call? */
    def isFinalCall = fqn contains "$final_call"

    def storeKey(workflow: OldStyleWorkflowDescriptor): ExecutionStoreKey = {
      finalCallCompanions.find(companion => fqn.endsWith("$final_call$" + companion.finalCallName)) match {
        case Some(companion) => FinalCallKey(companion.createCall(workflow))
        case None => throw new IllegalArgumentException(s"Cannot convert $fqn into a FinalCall store key for $workflow")
      }
    }
  }
}

/** Scope representing a "final call" that is inserted after all other workflow executions. */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait OldStyleFinalCall extends Scope {
  val companion: FinalCallCompanion[_ <: OldStyleFinalCall]
  val handle: FinalCallHandle
  val parent: Option[Scope] = Option(this.rootWorkflow)

  def workflow: OldStyleWorkflowDescriptor

  def execute(workflowMetadataResponse: WorkflowMetadataResponse)(implicit ec: ExecutionContext): Future[Unit]

  def poll(implicit ec: ExecutionContext, current: OldStyleExecutionHandle): Future[OldStyleExecutionHandle] = current match {
    case h if h == handle => Future.successful(handle)
    case _ => Future.failed(new IllegalStateException(s"Can only be polled with $handle."))
  }

  def prerequisiteCallNames: Set[LocallyQualifiedName] = {
    throw new UnsupportedOperationException("prerequisiteCallNames not supported for FinalCalls")
  }

  /**
    * This is handled as a special case by the WorkflowActor. We don't have to list the specially here.
    */
  def prerequisiteScopes: Set[Scope] = Set.empty

  override def rootWorkflow: Workflow = workflow.namespace.workflow

  override def unqualifiedName = "%s.$final_call$%s".format(workflow.name, companion.finalCallName)
}

/**
  * The companion object for each FinalCall allows instantiation of the class from the singleton, and also includes a
  * unique locally qualified name for each FinalCall type.
  *
  * @tparam A The type of the final call.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait FinalCallCompanion[A <: OldStyleFinalCall] {
  val finalCallName: String

  def createCall(workflow: OldStyleWorkflowDescriptor): A
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class FinalCallHandle extends OldStyleExecutionHandle {
  override val isDone = true
  override val result = OldStyleSuccessfulFinalCallExecution
}
