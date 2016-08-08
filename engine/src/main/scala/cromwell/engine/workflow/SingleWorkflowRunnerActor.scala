package cromwell.engine.workflow

import java.nio.file.Path
import java.util.UUID

import akka.actor.FSM.{CurrentState, Transition}
import akka.actor._
import akka.pattern.pipe
import akka.util.Timeout
import better.files._
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{WorkflowId, ExecutionStore => _, _}
import cromwell.engine.workflow.SingleWorkflowRunnerActor._
import cromwell.engine.workflow.WorkflowManagerActor.RetrieveNewWorkflows
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.SubmitWorkflow
import cromwell.server.CromwellRootActor
import cromwell.services.metadata.MetadataService.{GetSingleWorkflowMetadataAction, GetStatus, WorkflowOutputs}
import cromwell.util.PromiseActor._
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.metadata.MetadataBuilderActor
import spray.http.StatusCodes
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Try}

object SingleWorkflowRunnerActor {
  def props(source: WorkflowSourceFiles, metadataOutputFile: Option[Path]): Props = {
    Props(new SingleWorkflowRunnerActor(source, metadataOutputFile))
  }

  sealed trait RunnerMessage
  // The message to actually run the workflow is made explicit so the non-actor Main can `ask` this actor to do the
  // running and collect a result.
  case object RunWorkflow extends RunnerMessage
  private case object IssuePollRequest extends RunnerMessage
  private case object IssueReply extends RunnerMessage

  sealed trait RunnerState
  case object NotStarted extends RunnerState
  case object RunningWorkflow extends RunnerState
  case object RequestingOutputs extends RunnerState
  case object RequestingMetadata extends RunnerState
  case object Done extends RunnerState

  final case class RunnerData(replyTo: Option[ActorRef] = None,
                              terminalState: Option[WorkflowState] = None,
                              id: Option[WorkflowId] = None,
                              failures: Seq[Throwable] = Seq.empty) {

    def addFailure(message: String): RunnerData = addFailure(new Throwable(message))

    def addFailure(e: Throwable): RunnerData = this.copy(failures = e +: failures)
  }

  implicit class EnhancedJsObject(val jsObject: JsObject) extends AnyVal {
    def state: WorkflowState = WorkflowState.fromString(jsObject.fields("status").asInstanceOf[JsString].value)
  }

  private val Tag = "SingleWorkflowRunnerActor"
}

/**
 * Designed explicitly for the use case of the 'run' functionality in Main. This Actor will start a workflow,
 * print out the outputs when complete and then shut down the actor system. Note that multiple aspects of this
 * are sub-optimal for future use cases where one might want a single workflow being run.
 */
class SingleWorkflowRunnerActor(source: WorkflowSourceFiles, metadataOutputPath: Option[Path])
  extends CromwellRootActor with LoggingFSM[RunnerState, RunnerData] {

  import SingleWorkflowRunnerActor._
  private val backoff = SimpleExponentialBackoff(1 second, 1 minute, 1.2)
  private implicit val system = context.system
  private implicit val timeout = Timeout(5 seconds)

  startWith(NotStarted, RunnerData())

  private def requestMetadata: State = {
    val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor), s"MetadataRequest-Workflow-${stateData.id.get}")
    metadataBuilder ! GetSingleWorkflowMetadataAction(stateData.id.get, None, None)
    goto (RequestingMetadata)
  }

  private def schedulePollRequest(): Unit = {
    context.system.scheduler.scheduleOnce(backoff.backoffMillis.millis, self, IssuePollRequest)
  }

  private def requestStatus(): Unit = {
    // This requests status via the metadata service rather than instituting an FSM watch on the underlying workflow actor.
    // Cromwell's eventual consistency means it isn't safe to use an FSM transition to a terminal state as the signal for
    // when outputs or metadata have stabilized.
    val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor), s"StatusRequest-Workflow-${stateData.id.get}-request-${UUID.randomUUID()}")
    metadataBuilder.askNoTimeout(GetStatus(stateData.id.get)) pipeTo self
  }

  private def issueReply: State = {
    self ! IssueReply
    goto (Done)
  }

  when (NotStarted) {
    case Event(RunWorkflow, data) =>
      log.info(s"$Tag: Submitting workflow")
      workflowStoreActor ! SubmitWorkflow(source)
      goto (RunningWorkflow) using data.copy(replyTo = Option(sender()))
  }

  when (RunningWorkflow) {
    case Event(WorkflowStoreActor.WorkflowSubmittedToStore(id), data) =>
      log.info(s"$Tag: Workflow submitted UUID($id)")
      // Since we only have a single workflow, force the WorkflowManagerActor's hand in case the polling rate is long
      workflowManagerActor ! RetrieveNewWorkflows
      schedulePollRequest()
      stay() using data.copy(id = Option(id))
    case Event(IssuePollRequest, data) =>
      data.id match {
        case None => schedulePollRequest()
        case _ => requestStatus()
      }
      stay()
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), data) if !jsObject.state.isTerminal =>
      schedulePollRequest()
      stay()
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), data) if jsObject.state == WorkflowSucceeded =>
      val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor))
      metadataBuilder ! WorkflowOutputs(data.id.get)
      goto(RequestingOutputs) using data.copy(terminalState = Option(WorkflowSucceeded))
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), data) if jsObject.state == WorkflowFailed =>
      val updatedData = data.copy(terminalState = Option(WorkflowFailed)).addFailure(s"Workflow ${data.id.get} transitioned to state Failed")
      // If there's an output path specified then request metadata, otherwise issue a reply to the original sender.
      val nextState = if (metadataOutputPath.isDefined) requestMetadata else issueReply
      nextState using updatedData
  }

  when (RequestingOutputs) {
    case Event(RequestComplete((StatusCodes.OK, outputs: JsObject)), _) =>
      outputOutputs(outputs)
      if (metadataOutputPath.isDefined) requestMetadata else issueReply
  }

  when (RequestingMetadata) {
    case Event(RequestComplete((StatusCodes.OK, metadata: JsObject)), _) =>
      outputMetadata(metadata)
      issueReply
  }

  when (Done) {
    case Event(IssueReply, data) =>
      data.terminalState foreach { state => log.info(s"$Tag workflow finished with status '$state'.") }
      data.failures foreach { e => log.error(e, e.getMessage) }

      val message = data.terminalState collect { case WorkflowSucceeded => () } getOrElse Status.Failure(data.failures.head)
      data.replyTo foreach  { _ ! message }
      stay()
  }

  onTransition {
    case NotStarted -> RunningWorkflow => schedulePollRequest()
  }

  private def failAndFinish(e: Throwable): State = {
    log.error(e, s"$Tag received Failure message: ${e.getMessage}")
    issueReply using stateData.addFailure(e)
  }

  whenUnhandled {
    // Handle failures for all WorkflowManagerFailureResponses generically.
    case Event(r: WorkflowManagerFailureResponse, data) => failAndFinish(r.failure)
    case Event(Failure(e), data) => failAndFinish(e)
    case Event(Status.Failure(e), data) => failAndFinish(e)
    case Event(RequestComplete((_, snap)), _) => failAndFinish(new RuntimeException(s"Unexpected API completion message: $snap"))
    case Event((CurrentState(_, _) | Transition(_, _, _)), _) =>
      // ignore uninteresting current state and transition messages
      stay()
    case Event(m, _) =>
      log.warning(s"$Tag: received unexpected message: $m")
      stay()
  }

  /**
    * Outputs the outputs to stdout, and then requests the metadata.
    */
  private def outputOutputs(outputs: JsObject): Unit = {
    println(outputs.prettyPrint)
  }

  private def outputMetadata(metadata: JsObject): Try[Unit] = {
    Try {
      val path = metadataOutputPath.get
      if (path.isDirectory) {
        log.error("Specified metadata path is a directory, should be a file: " + path)
      } else {
        log.info(s"$Tag writing metadata to $path")
        path.createIfNotExists().write(metadata.prettyPrint)
      }
    }
  }
}
