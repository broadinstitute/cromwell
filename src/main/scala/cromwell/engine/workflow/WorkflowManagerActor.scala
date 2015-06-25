package cromwell.engine.workflow

import java.util.UUID

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.{ask, pipe}
import cromwell.binding
import cromwell.binding._
import cromwell.engine._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.workflow.WorkflowActor.{Restart, GetOutputs, Start}
import cromwell.util.WriteOnceStore
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.language.postfixOps


object WorkflowManagerActor {
  class WorkflowNotFoundException extends RuntimeException

  sealed trait WorkflowManagerActorMessage
  case class SubmitWorkflow(wdlSource: WdlSource, wdlJson: WdlJson, inputs: binding.WorkflowRawInputs) extends WorkflowManagerActorMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerActorMessage
  case object Shutdown extends WorkflowManagerActorMessage
  case class SubscribeToWorkflow(id: WorkflowId) extends WorkflowManagerActorMessage

  def props(dataAccess: DataAccess): Props = Props(new WorkflowManagerActor(dataAccess))

  case class RestartableWorkflow(id: WorkflowId, source: WdlSource, json: WdlJson, inputs: binding.WorkflowRawInputs)
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

  override def preStart() {
    restartIncompleteWorkflows()
  }

  def receive = LoggingReceive {
    case SubmitWorkflow(wdlSource, wdlJson, inputs) =>
      submitWorkflow(wdlSource, wdlJson, inputs, isRestart = false, maybeWorkflowId = None) pipeTo sender
    case WorkflowStatus(id) => dataAccess.getWorkflowState(id) pipeTo sender
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

  private def submitWorkflow(wdlSource: WdlSource, wdlJson: WdlJson, inputs: binding.WorkflowRawInputs,
                             isRestart: Boolean, maybeWorkflowId: Option[WorkflowId]): Future[WorkflowId] = {
    for {
      namespace <- Future(WdlNamespace.load(wdlSource))
      coercedInputs <- Future.fromTry(namespace.coerceRawInputs(inputs))
      descriptor = new WorkflowDescriptor(namespace, wdlSource, wdlJson, coercedInputs)
      workflowActor = context.actorOf(WorkflowActor.props(descriptor, backend, dataAccess))
      workflowId = maybeWorkflowId.getOrElse(UUID.randomUUID())
      _ <- Future.fromTry(workflowStore.insert(workflowId, workflowActor))
    } yield {
      workflowActor ! (if (isRestart) Restart else Start)
      workflowActor ! SubscribeTransitionCallBack(self)
      workflowId
    }
  }

  private def restartWorkflow(restartableWorkflow: RestartableWorkflow): Future[WorkflowId] = {
    submitWorkflow(restartableWorkflow.source, restartableWorkflow.json, restartableWorkflow.inputs,
      isRestart = true, Option(restartableWorkflow.id))
  }

  private def updateWorkflowState(workflow: WorkflowActorRef, state: WorkflowState): Future[Unit] = {
    val id = idByWorkflow(workflow)
    dataAccess.updateWorkflowState(id, state)
  }

  private def idByWorkflow(workflow: WorkflowActorRef): WorkflowId = {
    workflowStore.toMap.collectFirst { case (k, v) if v == workflow => k }.get
  }

  private def restartIncompleteWorkflows(): Unit = {
    // If the clob inputs for this workflow can be converted to JSON, return the JSON
    // version of those inputs in a Some().  Otherwise return None.
    def clobToJsonInputs(workflowInfo: WorkflowInfo): Option[binding.WorkflowRawInputs] = {
      workflowInfo.wdlJson.parseJson match {
        case JsObject(rawInputs) => Option(rawInputs)
        case x =>
          log.error(s"Error restarting workflow ${workflowInfo.workflowId}: expected JSON inputs, got '$x'.")
          None
      }
    }

    type QuantifierAndPlural = (String, String)
    def pluralize(num: Int): QuantifierAndPlural = {
      (if (num == 0) "no" else num.toString, if (num == 1) "" else "s")
    }

    val result = for {
      workflowInfos <- dataAccess.getWorkflowsByState(Seq(WorkflowSubmitted, WorkflowRunning))
    } yield {
        val restartableWorkflows = (for {
          workflowInfo <- workflowInfos
          jsonInputs = clobToJsonInputs(workflowInfo)
          if jsonInputs.isDefined
        } yield RestartableWorkflow(workflowInfo.workflowId, workflowInfo.wdlSource, workflowInfo.wdlJson, jsonInputs.get)).toSeq

        val num = restartableWorkflows.length
        val (displayNum, plural) = pluralize(num)
        log.info(s"Found $displayNum workflow$plural to restart.")

        if (num > 0) {
          val ids = restartableWorkflows.map { _.id.toString }.sorted
          log.info(s"Restarting workflow ID$plural: " + ids.mkString(", "))
        }

        restartableWorkflows foreach restartWorkflow
      }
    result recover {
      case e: Throwable => log.error(e, e.getMessage)
    }
  }
}
