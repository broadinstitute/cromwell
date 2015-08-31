package cromwell.engine.workflow

import akka.actor.FSM.{CurrentState, Transition}
import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import akka.pattern.ask
import cromwell.binding
import cromwell.binding.{WdlJson, WdlSource}
import cromwell.engine._
import cromwell.engine.workflow.WorkflowManagerActor.{SubmitWorkflow, SubscribeToWorkflow, WorkflowOutputs}
import spray.json._

import scala.concurrent.Await
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

object SingleWorkflowRunnerActor {
  def props(wdlSource: WdlSource, wdlJson: WdlJson, inputs: binding.WorkflowRawInputs, workflowManager: ActorRef): Props = {
    Props(new SingleWorkflowRunnerActor(wdlSource, wdlJson, inputs, workflowManager))
  }
}

/**
 * Designed explicitly for the use case of the 'run' functionality in Main. This Actor will start a workflow,
 * print out the outputs when complete and then shut down the actor system. Note that multiple aspects of this
 * are sub-optimal for future use cases where one might want a single workflow being run.
 */
case class SingleWorkflowRunnerActor(wdlSource: WdlSource,
                                     wdlJson: WdlJson,
                                     inputs: binding.WorkflowRawInputs,
                                     workflowManager: ActorRef) extends Actor with CromwellActor {
  val log = Logging(context.system, classOf[SingleWorkflowRunnerActor])
  val tag = "SingleWorkflowRunnerActor"
  // Note that id isn't used until *after* the submitWorkflow Future is complete
  private var id: WorkflowId = _

  override def preStart(): Unit = {
    log.info(s"$tag: launching workflow")
    val eventualId = workflowManager.ask(SubmitWorkflow(wdlSource, wdlJson, inputs)).mapTo[WorkflowId]
    eventualId onComplete {
      case Success(x) => subscribeToWorkflow(x)
      case Failure(e) =>
        log.error(e, s"$tag: ${e.getMessage}")
        terminate()
    }
  }

  def receive = {
    case Transition(_, _, state: WorkflowState) if state.isTerminal => handleTermination(state)
    case Transition(_, _, state: WorkflowState) => log.info(s"$tag: transitioning to $state")
    case m => log.warning(s"$tag: received unexpected message: $m")
  }

  private def handleTermination(state: WorkflowState): Unit = {
    log.info(s"$tag: workflow finished with status '$state'.")

    // If this is a successful termination, retrieve & print out the outputs
    if (state == WorkflowSucceeded) {
      val eventualOutputs = workflowManager.ask(WorkflowOutputs(id)).mapTo[binding.WorkflowOutputs]
      val outputs = Await.result(eventualOutputs, 5 seconds)
      import cromwell.binding.values.WdlValueJsonFormatter._
      println(outputs.toJson.prettyPrint)
    }

    terminate()
  }

  /**
   * Does a few things by side effect:
   *     - Sets 'id' to the submitted workflow id
   *     - Logs the workflow id
   *     - Subscribes to the new workflow's transition events
   */
  private def subscribeToWorkflow(workflowId: WorkflowId): Unit = {
    id = workflowId
    log.info(s"SingleWorkflowRunnerActor: workflow ID UUID($id)")
    workflowManager ! SubscribeToWorkflow(id)
  }

  private def terminate(): Unit = {
    context.system.shutdown()
  }
}
