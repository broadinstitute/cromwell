package cromwell

import akka.actor.ActorSystem
import cromwell.core.{CallOutput, WorkflowId}
import cromwell.engine.ExecutionStatus._
import cromwell.engine.backend.WorkflowDescriptor
import cromwell.engine.db.DataAccess.WorkflowExecutionAndAux
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.db.slick.Execution
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor.{MaterializationFailure, MaterializationSuccess, MaterializeWorkflow}
import org.joda.time.DateTime
import wdl4s._
import wdl4s.values.WdlValue

import scala.concurrent.{ExecutionContext, Future}
import scala.language.implicitConversions

package object engine {
  /**
   * Represents the collection of source files that a user submits to run a workflow
   */
  final case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  final case class AbortFunction(function: () => Unit)
  final case class AbortRegistrationFunction(register: AbortFunction => Unit)

  final case class ExecutionEventEntry(description: String, startTime: DateTime, endTime: DateTime)
  final case class QualifiedFailureEventEntry(workflowId: String, execution: Option[ExecutionDatabaseKey], failure: String, timestamp: DateTime) {
    def dequalify = FailureEventEntry(failure, timestamp)
  }
  final case class FailureEventEntry(failure: String, timestamp: DateTime)
  final case class CallAttempt(fqn: FullyQualifiedName, attempt: Int)

  type WorkflowOptionsJson = String
  type WorkflowOutputs = Map[FullyQualifiedName, CallOutput]
  type FullyQualifiedName = String

  type HostInputs = Map[String, WdlValue]

  class CromwellFatalException(exception: Throwable) extends Exception(exception)
  class PreemptedException(msg: String) extends Exception(msg)

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def isScatter = fqn.contains(Scatter.FQNIdentifier)
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, CallOutput]) extends AnyVal {
    def mapToValues: Map[A, WdlValue] = m map {
      case (k, CallOutput(wdlValue, hash)) => (k, wdlValue)
    }
  }

  implicit class EnhancedExecution(val execution: Execution) extends AnyVal {
    import cromwell.engine.ExecutionIndex._

    def isShard = execution.index.toIndex.isShard
    def isScatter = execution.callFqn.contains(Scatter.FQNIdentifier)
    def isCollector(keys: Traversable[Execution]): Boolean = {
      !execution.isShard &&
        (keys exists { e => (e.callFqn == execution.callFqn) && e.isShard })
    }
    def toKey: ExecutionDatabaseKey = ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex, execution.attempt)
    def executionStatus: ExecutionStatus = ExecutionStatus.withName(execution.status)
  }

  // Used to convert the database returned value `executionAndAux` to a WorkflowDescriptor
  def workflowDescriptorFromExecutionAndAux(executionAndAux: WorkflowExecutionAndAux)(implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[WorkflowDescriptor] = {
    //imports for implicits
    import cromwell.engine.db.slick.SlickDataAccess.ClobToRawString
    import cromwell.util.PromiseActor.EnhancedActorRef

    val id = WorkflowId.fromString(executionAndAux.execution.workflowExecutionUuid)
    val sources = WorkflowSourceFiles(executionAndAux.aux.wdlSource.toRawString, executionAndAux.aux.jsonInputs.toRawString, executionAndAux.aux.workflowOptions.toRawString)

    val materializeWorkflowDescriptorActor = actorSystem.actorOf(MaterializeWorkflowDescriptorActor.props())

    materializeWorkflowDescriptorActor.askNoTimeout(MaterializeWorkflow(id, sources))  map {
      case MaterializationSuccess(workflowDescriptor) => workflowDescriptor
      case MaterializationFailure(error) => throw error
    } andThen {
      case _ => actorSystem.stop(materializeWorkflowDescriptorActor)
    }
  }
}
