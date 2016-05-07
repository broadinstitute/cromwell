package cromwell

import akka.actor.ActorSystem
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.{CallOutput, WorkflowId}
import cromwell.engine.backend.OldStyleWorkflowDescriptor
import cromwell.engine.db.DataAccess.WorkflowExecutionAndAux
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.workflow.OldStyleMaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.OldStyleMaterializeWorkflowDescriptorActor.{MaterializeWorkflow, MaterializeWorkflowDescriptorFailure, MaterializeWorkflowDescriptorSuccess}
import org.joda.time.DateTime
import wdl4s._
import wdl4s.values.WdlValue

import scala.concurrent.{ExecutionContext, Future}
import scala.language.implicitConversions
import scala.util.{Failure, Success, Try}

package object engine {
  /**
   * Represents the collection of source files that a user submits to run a workflow
   */
  final case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  final case class AbortFunction(function: () => Unit)
  final case class AbortRegistrationFunction(register: AbortFunction => Unit)

  final case class QualifiedFailureEventEntry(workflowId: String, execution: Option[ExecutionDatabaseKey], failure: String, timestamp: DateTime) {
    def dequalify = FailureEventEntry(failure, timestamp)
  }
  final case class FailureEventEntry(failure: String, timestamp: DateTime)
  final case class CallAttempt(fqn: FullyQualifiedName, attempt: Int)

  type WorkflowOptionsJson = String
  type WorkflowOutputs = Map[FullyQualifiedName, CallOutput]
  type FullyQualifiedName = String

  type HostInputs = Map[String, WdlValue]

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
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

  object WorkflowFailureMode {
    def tryParse(mode: String): Try[WorkflowFailureMode] = {
      val modes = Seq(ContinueWhilePossible, NoNewCalls)
      modes find { _.toString.equalsIgnoreCase(mode) } map { Success(_) } getOrElse Failure(new Exception(s"Invalid workflow failure mode: $mode"))
    }
  }
  sealed trait WorkflowFailureMode {
    def allowNewCallsAfterFailure: Boolean
  }
  case object ContinueWhilePossible extends WorkflowFailureMode { override val allowNewCallsAfterFailure = true }
  case object NoNewCalls extends WorkflowFailureMode { override val allowNewCallsAfterFailure = false }

  // Used to convert the database returned value `executionAndAux` to a WorkflowDescriptor
  def workflowDescriptorFromExecutionAndAux(executionAndAux: WorkflowExecutionAndAux)(implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[OldStyleWorkflowDescriptor] = {
    //imports for implicits
    import cromwell.database.SqlConverters.ClobToRawString
    import cromwell.util.PromiseActor.EnhancedActorRef

    val id = WorkflowId.fromString(executionAndAux.execution.workflowExecutionUuid)
    val sources = WorkflowSourceFiles(executionAndAux.aux.wdlSource.toRawString, executionAndAux.aux.jsonInputs.toRawString, executionAndAux.aux.workflowOptions.toRawString)

    val materializeWorkflowDescriptorActor = actorSystem.actorOf(OldStyleMaterializeWorkflowDescriptorActor.props())

    materializeWorkflowDescriptorActor.askNoTimeout(MaterializeWorkflow(id, sources))  map {
      case MaterializeWorkflowDescriptorSuccess(workflowDescriptor) => workflowDescriptor
      case MaterializeWorkflowDescriptorFailure(error) => throw error
    } andThen {
      case _ => actorSystem.stop(materializeWorkflowDescriptorActor)
    }
  }

  final case class EngineWorkflowDescriptor(backendDescriptor: BackendWorkflowDescriptor,
                                            declarations: WorkflowCoercedInputs,
                                            backendAssignments: Map[Call, String],
                                            failureMode: WorkflowFailureMode) {
    def id = backendDescriptor.id
    def namespace = backendDescriptor.workflowNamespace
  }
}
