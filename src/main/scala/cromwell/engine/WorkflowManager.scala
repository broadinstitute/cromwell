package cromwell.engine

import java.util.UUID

import akka.actor.{Props, Actor, ActorSystem, ActorRef}
import akka.pattern.pipe
import cromwell.util.WriteOnceStore
import scala.concurrent.ExecutionContext.Implicits.global
import cromwell.binding.{WdlBinding, WorkflowInputs, WdlSource}
import cromwell.engine.WorkflowActor.Start
import cromwell.engine.WorkflowManagerActor.{SubmitWorkflow, WorkflowStatus}
import cromwell.engine.backend.Backend
import cromwell.engine.backend.local.LocalBackend

import scala.concurrent.Future
import scala.util.Try

/**
 * Abstract notion of a workflow manager. Mainly exists to allow for actor-system-free unit testing of basic concepts
 */
trait WorkflowManager {
  type Workflow

  case class ManagedWorkflow(id: WorkflowId, workflow: Workflow)

  val backend: Backend

  def generateWorkflow(id: WorkflowId, wdl: WdlSource, inputs: WorkflowInputs): Try[Workflow]
  def workflowStatus(id: WorkflowId): Option[WorkflowState]

  private val workflowStore = new WriteOnceStore[WorkflowId, Workflow]

  /**
   * Generates a workflow with an ID, inserts it into the store and returns the pair.
   */
  def submitWorkflow(wdl: WdlSource, inputs: WorkflowInputs): Try[ManagedWorkflow] = {
    val workflowId = UUID.randomUUID()
    for {
      wf <- generateWorkflow(workflowId, wdl, inputs)
      x <- workflowStore.insert(workflowId, wf) // Come for the side effect, stay for the Try
    } yield ManagedWorkflow(workflowId, wf)
  }

  def workflowById(id: WorkflowId): Option[Workflow] = workflowStore.toMap.get(id)
}

/**
 * A WorkflowManager using WorkflowActors as the workflow handler
 */
trait ActorWorkflowManager extends WorkflowManager {
  override type Workflow = ActorRef // TODO: In a world where Akka Typed is no longer experimental switch to that
  override val backend = new LocalBackend
  implicit def actorSystem: ActorSystem

  override def generateWorkflow(id: WorkflowId, wdl: WdlSource, inputs: WorkflowInputs): Try[Workflow] = {
    val binding = Try(WdlBinding.process(wdl))
    binding map {x => actorSystem.actorOf(WorkflowActor.buildWorkflowActorProps(id, x, inputs))}
  }

  override def workflowStatus(id: WorkflowId): Option[WorkflowState] = workflowById(id) map {x => workflowState(x)}

  override def submitWorkflow(wdl: WdlSource, inputs: WorkflowInputs): Try[ManagedWorkflow] = {
    val managedWorkflow = super.submitWorkflow(wdl, inputs)
    managedWorkflow foreach {_.workflow ! Start(backend)}
    managedWorkflow
  }

  // FIXME: Placeholder for now because there's no facility in a WorkflowActor to get a state
  private def workflowState(workflow: ActorRef): WorkflowState = WorkflowRunning
}

object WorkflowManagerActor {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage

  def props(actorSystem: ActorSystem): Props = Props(classOf[WorkflowManagerActor], actorSystem)
}

/**
 * Responses to messages:
 * SubmitWorkflow: Returns a Future[Try[ManagedWorkflow]]
 * WorkflowStatus: Returns a Future[Option[WorkflowState]]
 *
 * FIXME: It might be nice to have a notion of a actor-based workflow manager which didn't have actor-based workflows but I doubt we'll care any time soon
 */
case class WorkflowManagerActor(actorSystem: ActorSystem) extends Actor with ActorWorkflowManager {
  def receive = {
    case SubmitWorkflow(wdl, inputs) => handleRequest(sender(), Future {submitWorkflow(wdl, inputs)})
    case WorkflowStatus(id) => handleRequest(sender(), Future {workflowStatus(id)})
  }

  private def handleRequest[T](sender: ActorRef, response: Future[T]): Unit = response pipeTo sender
}
