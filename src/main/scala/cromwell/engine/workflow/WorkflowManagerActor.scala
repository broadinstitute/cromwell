package cromwell.engine.workflow

import java.util.UUID

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.{ask, pipe}
import cromwell.binding
import cromwell.binding.{WdlNamespace, WdlSource, WorkflowDescriptor}
import cromwell.engine._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.{DataAccess, QueryWorkflowExecutionResult}
import cromwell.engine.workflow.WorkflowActor.{Start, GetOutputs}
import cromwell.util.WriteOnceStore
import spray.json._

import scala.collection.concurrent.TrieMap
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.language.postfixOps


object WorkflowManagerActor {
  class WorkflowNotFoundException extends RuntimeException

  sealed trait ActorWorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs) extends ActorWorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends ActorWorkflowManagerMessage
  case class WorkflowOutputs(id: WorkflowId) extends ActorWorkflowManagerMessage
  case object Shutdown extends ActorWorkflowManagerMessage
  case class SubscribeToWorkflow(id: WorkflowId) extends ActorWorkflowManagerMessage

  def props(dataAccess: DataAccess): Props = Props(new WorkflowManagerActor(dataAccess))
}

/**
 * Responses to messages:
 * SubmitWorkflow: Returns a `Future[WorkflowId]`
 * WorkflowStatus: Returns a `Future[Option[WorkflowState]]`
 * WorkflowOutputs: Returns a `Future[Option[binding.WorkflowOutputs]]` aka `Future[Option[Map[String, WdlValue]]`
 *
 */
class WorkflowManagerActor(dataAccess: DataAccess) extends Actor with CromwellActor {
  import WorkflowManagerActor._
  private val log = Logging(context.system, this)

  type WorkflowActorRef = ActorRef

  private val backend = new LocalBackend
  private val workflowStore = new WriteOnceStore[WorkflowId, WorkflowActorRef]
  // This *should* be persisted
  private val workflowStates = TrieMap.empty[WorkflowId, WorkflowState]

  override def preStart() {
    restartIncompleteWorkflows()
  }

  def receive = LoggingReceive {
    case SubmitWorkflow(wdl, inputs) => submitWorkflow(wdl, inputs) pipeTo sender
    case WorkflowStatus(id) => sender ! workflowStates.get(id)
    case Shutdown => context.system.shutdown()
    case WorkflowOutputs(id) => workflowOutputs(id) pipeTo sender
    case CurrentState(actor, state: WorkflowState) => updateWorkflowState(actor, state)
    case Transition(actor, oldState, newState: WorkflowState) => updateWorkflowState(actor, newState)
    case SubscribeToWorkflow(id) =>
      //  NOTE: This fails silently. Currently we're ok w/ this, but you might not be in the future
      workflowStore.toMap.get(id) foreach {_ ! SubscribeTransitionCallBack(sender())}
  }

  private def workflowOutputs(id: WorkflowId): Future[binding.WorkflowOutputs] = {
   workflowStore.toMap.get(id) map workflowToOutputs getOrElse Future.failed(new WorkflowNotFoundException)
  }

  private def workflowToOutputs(workflow: WorkflowActorRef): Future[binding.WorkflowOutputs] = {
    workflow.ask(GetOutputs).mapTo[binding.WorkflowOutputs]
  }

  /**
   * Processes input WDL, creates and starts a workflow actor, inserts it into the store and returns its ID.
   */
  def submitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs, workflowId: UUID = UUID.randomUUID()): Future[WorkflowId] = {
    for {
      namespace <- Future(WdlNamespace.load(wdl))
      coercedInputs <- Future.fromTry(namespace.coerceRawInputs(inputs))
      descriptor = new WorkflowDescriptor(namespace, coercedInputs)
      workflowActor = context.actorOf(WorkflowActor.props(descriptor, backend))
      _ <- Future.fromTry(workflowStore.insert(workflowId, workflowActor))
    } yield {
      workflowStates.put(workflowId, WorkflowSubmitted)
      workflowActor ! Start
      workflowActor ! SubscribeTransitionCallBack(self)
      workflowId
    }
  }

  private def updateWorkflowState(workflow: WorkflowActorRef, state: WorkflowState): Unit = {
    idByWorkflow(workflow) foreach { id => workflowStates.put(id, state) }
  }

  private def idByWorkflow(workflow: WorkflowActorRef): Option[WorkflowId] = {
    workflowStore.toMap collectFirst {case (k, v) if v == workflow => k}
  }

  private def restartIncompleteWorkflows(): Unit = {
    // If the clob inputs for this workflow can be converted to JSON, return the JSON
    // version of those inputs in a Some().  Otherwise return None.
    def clobToJsonInputs(result: QueryWorkflowExecutionResult): Option[binding.WorkflowRawInputs] = {
      result.wdlRawInputs.parseJson match {
        case JsObject(rawInputs) => Option(rawInputs)
        case x =>
          log.error(s"Error restarting workflow ${result.workflowId}: expected JSON inputs, got '$x'")
          None
      }
    }
    val RestartableStates = Some(Seq(WorkflowSubmitted, WorkflowRunning))
    // Attempt to restart all the workflows in restartable states whose clob raw inputs
    // can successfully be converted to JSON.
    val restartableWorkflows = for {
      workflow <- dataAccess.query(states = RestartableStates)
      jsonInputs = clobToJsonInputs(workflow)
      if jsonInputs.isDefined
    } yield (workflow, jsonInputs.get)

    restartableWorkflows foreach { case (workflow, jsonInputs) =>
      submitWorkflow(workflow.wdlSource, jsonInputs, workflow.workflowId)
    }
  }
}
