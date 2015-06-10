package cromwell.engine

import java.util.UUID

import akka.actor.FSM.{Transition, CurrentState, SubscribeTransitionCallBack}
import akka.actor.{Actor, ActorRef, Props}
import akka.event.LoggingReceive
import akka.pattern.{ask, pipe}
import cromwell.binding
import cromwell.binding.{WdlNamespace, WdlSource}
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.WriteOnceStore

import scala.collection.concurrent.TrieMap
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future


object WorkflowManagerActor {
  class WorkflowNotFoundException extends RuntimeException

  sealed trait ActorWorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs) extends ActorWorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends ActorWorkflowManagerMessage
  case class WorkflowOutputs(id: WorkflowId) extends ActorWorkflowManagerMessage
  case object Shutdown extends ActorWorkflowManagerMessage
  case class SubscribeToWorkflow(id: WorkflowId) extends ActorWorkflowManagerMessage

  def props: Props = Props(classOf[WorkflowManagerActor])
}

/**
 * Responses to messages:
 * SubmitWorkflow: Returns a `Future[WorkflowId]`
 * WorkflowStatus: Returns a `Future[Option[WorkflowState]]`
 * WorkflowOutputs: Returns a `Future[Option[binding.WorkflowOutputs]]` aka `Future[Option[Map[String, WdlValue]]`
 *
 */
class WorkflowManagerActor extends Actor {
  import WorkflowManagerActor._

  type WorkflowActorRef = ActorRef

  private val backend = new LocalBackend
  private val workflowStore = new WriteOnceStore[WorkflowId, WorkflowActorRef]
  // This *should* be persisted
  private val workflowStates = TrieMap.empty[WorkflowId, WorkflowState]

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
  private def submitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs): Future[WorkflowId] = {
    val workflowId = UUID.randomUUID()
    for {
      namespace <- Future(WdlNamespace.load(wdl))
      coercedInputs <- Future.fromTry(namespace.coerceRawInputs(inputs))
      workflowActor = context.actorOf(WorkflowActor.props(workflowId, namespace, coercedInputs, backend))
      _ <- Future.fromTry(workflowStore.insert(workflowId, workflowActor))
    } yield {
      workflowStates.put(workflowId, WorkflowSubmitted)
      workflowActor ! Start
      workflowActor ! SubscribeTransitionCallBack(self)
      workflowId
    }
  }

  private def updateWorkflowState(workflow: WorkflowActorRef, state: WorkflowState): Unit = {
    idByWorkflow(workflow) map { w => workflowStates.put(w, state) }
  }

  private def idByWorkflow(workflow: WorkflowActorRef): Option[WorkflowId] = {
    workflowStore.toMap collectFirst {case (k, v) if v == workflow => k}
  }
}
