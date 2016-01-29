package cromwell.engine.workflow

import java.nio.file.Path

import akka.actor.FSM.Transition
import akka.actor._
import better.files._
import cromwell.engine
import cromwell.engine.{CallOutput, _}
import cromwell.engine.workflow.SingleWorkflowRunnerActor._
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.{WdlValueJsonFormatter, WorkflowMetadataResponse}
import spray.json._

import scala.util._

object SingleWorkflowRunnerActor {
  def props(source: WorkflowSourceFiles, metadataOutputFile: Option[Path], workflowManager: ActorRef): Props = {
    Props(classOf[SingleWorkflowRunnerActor], source, metadataOutputFile, workflowManager)
  }

  sealed trait RunnerMessage
  // The message to actually run the workflow is made explicit so the non-actor Main can `ask` this actor to do the
  // running and collect a result.
  case object RunWorkflow extends RunnerMessage
  private case object IssueReply extends RunnerMessage

  sealed trait RunnerState
  case object NotStarted extends RunnerState
  case object RunningWorkflow extends RunnerState
  case object RequestingOutputs extends RunnerState
  case object RequestingMetadata extends RunnerState
  case object Done extends RunnerState

  final case class RunnerData(replyTo: Option[ActorRef] = None,
                              id: Option[WorkflowId] = None,
                              terminalState: Option[WorkflowState] = None,
                              failures: Seq[Throwable] = Seq.empty) {

    def addFailure(message: String): RunnerData = addFailure(new Throwable(message))

    def addFailure(e: Throwable): RunnerData = this.copy(failures = e +: failures)
  }
}

/**
 * Designed explicitly for the use case of the 'run' functionality in Main. This Actor will start a workflow,
 * print out the outputs when complete and then shut down the actor system. Note that multiple aspects of this
 * are sub-optimal for future use cases where one might want a single workflow being run.
 */
case class SingleWorkflowRunnerActor(source: WorkflowSourceFiles,
                                     metadataOutputPath: Option[Path],
                                     workflowManager: ActorRef) extends LoggingFSM[RunnerState, RunnerData] with CromwellActor {

  val tag = "SingleWorkflowRunnerActor"

  startWith(NotStarted, RunnerData())

  private def requestMetadata: State = {
    workflowManager ! WorkflowMetadata(stateData.id.get)
    goto (RequestingMetadata)
  }

  private def issueReply: State = {
    self ! IssueReply
    goto (Done)
  }

  when (NotStarted) {
    case Event(RunWorkflow, data) =>
      log.info(s"$tag: launching workflow")
      workflowManager ! SubmitWorkflow(source)
      goto (RunningWorkflow) using data.copy(replyTo = Option(sender()))
  }

  when (RunningWorkflow) {
    case Event(WorkflowManagerSubmitSuccess(id), data) =>
      log.info(s"$tag: workflow ID UUID($id)")
      workflowManager ! SubscribeToWorkflow(id)

      if (getAbortJobsOnTerminate) {
        Runtime.getRuntime.addShutdownHook(new Thread() {
          override def run(): Unit = {
            workflowManager ! WorkflowAbort(id)
            log.info(s"$tag: Waiting for workflow $id to abort...")
            while(stateName != Done)
              Thread.sleep(1000)
            log.info(s"$tag: Workflow $id aborted.")
          }
        })
      }
      stay using data.copy(id = Option(id))
    case Event(Transition(_, _, WorkflowSucceeded), data) =>
      workflowManager ! WorkflowOutputs(data.id.get)
      goto(RequestingOutputs) using data.copy(terminalState = Option(WorkflowSucceeded))
    case Event(Transition(_, _, state: WorkflowState), data) if state.isTerminal =>
      // A terminal state that is not `WorkflowSucceeded` is a failure.
      val updatedData = data.copy(terminalState = Option(state)).addFailure(s"Workflow ${data.id.get} transitioned to state $state")
      // If there's an output path specified then request metadata, otherwise issue a reply to the original sender.
      val nextState = if (metadataOutputPath.isDefined) requestMetadata else issueReply
      nextState using updatedData
    case Event(Transition(_, _, state), _) =>
      log.info(s"$tag: transitioning to $state")
      stay()
  }

  when (RequestingOutputs) {
    case Event(WorkflowManagerWorkflowOutputsSuccess(id, outputs), data) =>
      // Outputs go to stdout
      outputOutputs(outputs)
      if (metadataOutputPath.isDefined) requestMetadata else issueReply
  }
  
  when (RequestingMetadata) {
    case Event(r: WorkflowManagerWorkflowMetadataSuccess, data) =>
      val updatedData = outputMetadata(r.response) match {
        case Success(_) => data
        case Failure(e) => data.addFailure(e)
      }
      issueReply using updatedData
  }

  when (Done) {
    case Event(IssueReply, data) =>
      data.terminalState foreach { state => log.info(s"$tag workflow finished with status '$state'.") }
      data.failures foreach { e => log.error(e, e.getMessage) }

      val message = data.terminalState collect { case WorkflowSucceeded => () } getOrElse Status.Failure(data.failures.head)
      data.replyTo foreach  { _ ! message }
      stay()
  }

  private def failAndFinish(e: Throwable): State = {
    log.error(e, s"$tag received Failure message: ${e.getMessage}")
    issueReply using stateData.addFailure(e)
  }

  whenUnhandled {
    // Handle failures for all WorkflowManagerFailureResponses generically.
    case Event(r: WorkflowManagerFailureResponse, data) => failAndFinish(r.failure)
    case Event(Failure(e), data) => failAndFinish(e)
    case Event(Status.Failure(e), data) => failAndFinish(e)
    case Event(m, _) =>
      log.warning(s"$tag: received unexpected message: $m")
      stay()
  }

  /**
    * Outputs the outputs to stdout, and then requests the metadata.
    */
  private def outputOutputs(outputs: engine.WorkflowOutputs): Unit = {
    import WdlValueJsonFormatter._
    val outputValues = outputs mapValues { case CallOutput(wdlValue, _) => wdlValue }
    println(outputValues.toJson.prettyPrint)
  }

  private def outputMetadata(metadata: WorkflowMetadataResponse): Try[Unit] = {
    Try {
      // This import is required despite what IntelliJ thinks.
      import cromwell.webservice.WorkflowJsonSupport._
      val path = metadataOutputPath.get
      log.info(s"$tag writing metadata to $path")
      path.createIfNotExists().write(metadata.toJson.prettyPrint)
    }
  }
}
