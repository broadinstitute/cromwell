package cromwell.engine

import java.util.UUID

import akka.actor.{Actor, ActorRef, Props}
import akka.event.LoggingReceive
import akka.pattern.{ask, pipe}
import cromwell.binding
import cromwell.binding.{WdlBinding, WdlSource}
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.WriteOnceStore

import scala.collection.concurrent.TrieMap
import scala.collection.mutable.{Set => MutableSet}
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.{Future, Promise}


object WorkflowManagerActor {
  sealed trait ActorWorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs) extends ActorWorkflowManagerMessage
  case class NotifyCompletion(caller: Promise[Unit]) extends ActorWorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends ActorWorkflowManagerMessage
  case class WorkflowOutputs(id: WorkflowId) extends ActorWorkflowManagerMessage
  case object Shutdown extends ActorWorkflowManagerMessage
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
  private val listeners = MutableSet.empty[Promise[Unit]]

  def receive = LoggingReceive {
    case SubmitWorkflow(wdl, inputs) =>
      submitWorkflow(wdl, inputs) pipeTo sender

    case WorkflowStatus(id) =>
      sender ! workflowStates.get(id)

    case NotifyCompletion(promise) =>
      listeners += promise

    case WorkflowOutputs(id) =>
      // FIXME: What if the workflow isn't done? How best to handle?
      workflowOutputs(id) pipeTo sender

    case Shutdown =>
      context.system.shutdown()

    case WorkflowActor.Started =>
      updateWorkflowState(sender(), WorkflowRunning)

    case WorkflowActor.Done(outputs) =>
      updateTerminalWorkflowState(WorkflowSucceeded)

    case WorkflowActor.Failed(failures) =>
      updateTerminalWorkflowState(WorkflowFailed)
  }

  private def updateTerminalWorkflowState(state: WorkflowState): Unit = {
    listeners foreach { listener =>
      if (!listener.isCompleted) listener.success(())
    }
    updateWorkflowState(sender(), state)
  }

  private def workflowOutputs(id: WorkflowId): Future[Option[binding.WorkflowOutputs]] = {
    workflowStore.toMap.get(id) map workflowToOutputs getOrElse Future{None}
  }

  private def workflowToOutputs(workflow: WorkflowActorRef): Future[Option[binding.WorkflowOutputs]] =
    (workflow ? GetOutputs).mapTo[binding.WorkflowOutputs] map { Option(_) }

  /**
   * Processes input WDL, creates and starts a workflow actor, inserts it into the store and returns its ID.
   */
  private def submitWorkflow(wdl: WdlSource, inputs: binding.WorkflowRawInputs): Future[WorkflowId] = {
    val workflowId = UUID.randomUUID()
    for {
      binding <- Future(WdlBinding.process(wdl))
      coercedInputs <- Future.fromTry(binding.coerceRawInputs(inputs))
      workflowActor = context.actorOf(WorkflowActor.props(workflowId, binding, coercedInputs, backend))
      _ <- Future.fromTry(workflowStore.insert(workflowId, workflowActor))
    } yield {
      workflowStates.put(workflowId, WorkflowSubmitted)
      workflowActor ! Start
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
