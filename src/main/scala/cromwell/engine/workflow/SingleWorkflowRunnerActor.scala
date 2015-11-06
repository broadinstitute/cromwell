package cromwell.engine.workflow

import java.nio.file.Path

import akka.actor.FSM.{CurrentState, Transition}
import akka.actor.{Actor, ActorRef, Props, Status}
import akka.event.Logging
import better.files._
import cromwell.binding
import cromwell.binding.{CallOutput, FullyQualifiedName}
import cromwell.binding.values.WdlValue
import cromwell.engine._
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.webservice.WorkflowMetadataResponse
import spray.json._

object SingleWorkflowRunnerActor {
  def props(source: WorkflowSourceFiles, metadataOutputFile: Option[Path], workflowManager: ActorRef): Props = {
    Props(classOf[SingleWorkflowRunnerActor], source, metadataOutputFile, workflowManager)
  }
}

/**
 * Designed explicitly for the use case of the 'run' functionality in Main. This Actor will start a workflow,
 * print out the outputs when complete and then shut down the actor system. Note that multiple aspects of this
 * are sub-optimal for future use cases where one might want a single workflow being run.
 */
case class SingleWorkflowRunnerActor(source: WorkflowSourceFiles,
                                     metadataOutputPath: Option[Path],
                                     workflowManager: ActorRef) extends Actor with CromwellActor {
  val log = Logging(context.system, classOf[SingleWorkflowRunnerActor])
  val tag = "SingleWorkflowRunnerActor"
  // Note that id isn't used until *after* the submitWorkflow Future is complete
  private var id: WorkflowId = _

  override def preStart(): Unit = {
    log.info(s"$tag: launching workflow")
    workflowManager ! SubmitWorkflow(source)
  }

  def receive = {
    case workflowId: WorkflowId => subscribeToWorkflow(workflowId)
    case outputs: Map[FullyQualifiedName@unchecked, CallOutput@unchecked] => outputOutputs(outputs)
    case metadata: WorkflowMetadataResponse => outputMetadata(metadata) // pkg "webservice", but this is what WMA sends
    case currentState: CurrentState[_] => log.debug(s"$tag: ignoring current state message: $currentState")
    case Transition(_, _, state: WorkflowState) if state.isTerminal => handleTermination(state)
    case Transition(_, _, state: WorkflowState) => log.info(s"$tag: transitioning to $state")
    case Status.Failure(t) => handleFailure(t) // .pipeTo() transforms util.Failure() to Status.Failure()
    case m => log.warning(s"$tag: received unexpected message: $m")
  }

  /**
   * Logs the workflow completion, then either a) if the workflow completed successfully: requests the outputs, and
   * later the metadata, or b) if not did not complete successfully, immediately requests the metadata.
   * @param state The workflow termination state.
   */
  private def handleTermination(state: WorkflowState): Unit = {
    log.info(s"$tag: workflow finished with status '$state'.")

    // If this is a successful termination, retrieve & print out the outputs
    if (state == WorkflowSucceeded) {
      workflowManager ! WorkflowOutputs(id)
    } else {
      sendMetadataMessage()
    }
  }

  /**
   * Logs the errors and shuts down.
   * @param t The error.
   */
  private def handleFailure(t: Throwable): Unit = {
    log.error(t, s"$tag: ${t.getMessage}")
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

  /**
   * Outputs the outputs to stdout, and then requests the metadata.
   */
  private def outputOutputs(outputs: binding.WorkflowOutputs): Unit = {
    import cromwell.binding.values.WdlValueJsonFormatter._
    val outputValues = outputs map {
      case (k, CallOutput(wdlValue, hash)) => (k, wdlValue)
    }
    println(outputValues.toJson.prettyPrint)
    sendMetadataMessage()
  }

  /**
   * Requests the metadata, if the path is specified, or shuts down.
   */
  private def sendMetadataMessage(): Unit = {
    metadataOutputPath match {
      case Some(file) => workflowManager ! WorkflowMetadata(id)
      case None => terminate()
    }
  }

  /**
   * Outputs the metadata, then shuts down.
   */
  private def outputMetadata(metadata: WorkflowMetadataResponse): Unit = {
    try {
      import cromwell.webservice.WorkflowJsonSupport._
      val path = metadataOutputPath.get
      log.info(s"Writing metadata to $path")
      path.createIfNotExists().write(metadata.toJson.prettyPrint)
      terminate()
    } catch {
      case e: Exception =>
        handleFailure(e)
    }
  }

  /**
   * Got to tell the manager its job is done.
   */
  private def terminate(): Unit = {
    workflowManager ! Shutdown
  }
}
