package cromwell.engine.workflow

import java.nio.file.Path

import akka.actor.FSM.{CurrentState, Transition}
import akka.actor._
import better.files._
import cromwell.core.{JobOutput, WorkflowId, _}
import cromwell.engine
import cromwell.engine._
import cromwell.engine.workflow.OldStyleSingleWorkflowRunnerActor._
import cromwell.engine.workflow.OldStyleWorkflowManagerActor.{WorkflowMetadata, _}
import cromwell.engine.workflow.WorkflowActor.WorkflowSucceededState
import cromwell.engine.workflow.WorkflowManagerActor.{SubmitWorkflowCommand, SubscribeToWorkflowCommand}
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.{WdlValueJsonFormatter, WorkflowMetadataResponse}
import spray.json._

import scala.util._
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleSingleWorkflowRunnerActor {
  def props(source: WorkflowSourceFiles, metadataOutputFile: Option[Path], workflowManager: ActorRef): Props = {
    Props(classOf[OldStyleSingleWorkflowRunnerActor], source, metadataOutputFile, workflowManager)
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
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleSingleWorkflowRunnerActor(source: WorkflowSourceFiles,
                                             metadataOutputPath: Option[Path],
                                             workflowManager: ActorRef) extends LoggingFSM[RunnerState, RunnerData] with CromwellActor {

  import OldStyleSingleWorkflowRunnerActor._

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
      log.info(s"$tag: Launching workflow")
      val submitMessage = SubmitWorkflowCommand(source)
      workflowManager ! submitMessage
      goto (RunningWorkflow) using data.copy(replyTo = Option(sender()))
  }

  when (RunningWorkflow) {
    case Event(WorkflowManagerSubmitSuccess(id), data) =>
      log.info(s"$tag: Workflow submitted UUID($id)")
      workflowManager ! SubscribeToWorkflowCommand(id)
      stay() using data.copy(id = Option(id))
    case Event(Transition(_, _, WorkflowSucceededState), data) =>
      workflowManager ! WorkflowOutputs(data.id.get)
      goto(RequestingOutputs) using data.copy(terminalState = Option(WorkflowSucceeded))
    case Event(Transition(_, _, WorkflowFailed), data) =>
      val updatedData = data.copy(terminalState = Option(WorkflowFailed)).addFailure(s"Workflow ${data.id.get} transitioned to state Failed")
      // If there's an output path specified then request metadata, otherwise issue a reply to the original sender.
      val nextState = if (metadataOutputPath.isDefined) requestMetadata else issueReply
      nextState using updatedData
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
    case Event((CurrentState(_, _) | Transition(_, _, _)), _) =>
      // ignore uninteresting current state and transition messages
      stay()
    case Event(m, _) =>
      log.warning(s"$tag: received unexpected message: $m")
      stay()
  }

  /**
    * Outputs the outputs to stdout, and then requests the metadata.
    */
  private def outputOutputs(outputs: engine.WorkflowOutputs): Unit = {
    import WdlValueJsonFormatter._
    val outputValues = outputs mapValues { case JobOutput(wdlValue, _) => wdlValue }
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
