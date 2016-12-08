package cromwell.engine.workflow

import java.nio.file.Path
import java.util.UUID

import akka.actor.FSM.{CurrentState, Transition}
import akka.actor._
import better.files._
import cats.instances.try_._
import cats.syntax.functor._
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.engine.workflow.SingleWorkflowRunnerActor._
import cromwell.engine.workflow.WorkflowManagerActor.RetrieveNewWorkflows
import cromwell.engine.workflow.workflowstore.{InMemoryWorkflowStore, WorkflowStoreActor}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.SubmitWorkflow
import cromwell.jobstore.EmptyJobStoreActor
import cromwell.server.CromwellRootActor
import cromwell.services.metadata.MetadataService.{GetSingleWorkflowMetadataAction, GetStatus, WorkflowOutputs}
import cromwell.subworkflowstore.EmptySubWorkflowStoreActor
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.metadata.MetadataBuilderActor
import spray.http.StatusCodes
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Try}

/**
 * Designed explicitly for the use case of the 'run' functionality in Main. This Actor will start a workflow,
 * print out the outputs when complete and reply with a result.
 */
class SingleWorkflowRunnerActor(source: WorkflowSourceFilesCollection, metadataOutputPath: Option[Path])
  extends CromwellRootActor with LoggingFSM[RunnerState, SwraData] {

  override val serverMode = false

  import SingleWorkflowRunnerActor._
  private val backoff = SimpleExponentialBackoff(1 second, 1 minute, 1.2)

  override val abortJobsOnTerminate = true
  override lazy val workflowStore = new InMemoryWorkflowStore()
  override lazy val jobStoreActor = context.actorOf(EmptyJobStoreActor.props)
  override lazy val subWorkflowStoreActor = context.actorOf(EmptySubWorkflowStoreActor.props)

  startWith(NotStarted, EmptySwraData)

  when (NotStarted) {
    case Event(RunWorkflow, EmptySwraData) =>
      log.info(s"$Tag: Submitting workflow")
      workflowStoreActor ! SubmitWorkflow(source)
      goto(SubmittedWorkflow) using SubmittedSwraData(sender())
  }

  when (SubmittedWorkflow) {
    case Event(WorkflowStoreActor.WorkflowSubmittedToStore(id), SubmittedSwraData(replyTo)) =>
      log.info(s"$Tag: Workflow submitted UUID($id)")
      // Since we only have a single workflow, force the WorkflowManagerActor's hand in case the polling rate is long
      workflowManagerActor ! RetrieveNewWorkflows
      schedulePollRequest()
      goto(RunningWorkflow) using RunningSwraData(replyTo, id)
  }

  when (RunningWorkflow) {
    case Event(IssuePollRequest, RunningSwraData(_, id)) =>
      requestStatus(id)
      stay()
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), RunningSwraData(_, _)) if !jsObject.state.isTerminal =>
      schedulePollRequest()
      stay()
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), RunningSwraData(replyTo, id)) if jsObject.state == WorkflowSucceeded =>
      val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor),
        s"CompleteRequest-Workflow-$id-request-${UUID.randomUUID()}")
      metadataBuilder ! WorkflowOutputs(id)
      log.info(s"$Tag workflow finished with status '$WorkflowSucceeded'.")
      goto(RequestingOutputs) using SucceededSwraData(replyTo, id)
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), RunningSwraData(replyTo, id)) if jsObject.state == WorkflowFailed =>
      log.info(s"$Tag workflow finished with status '$WorkflowFailed'.")
      requestMetadataOrIssueReply(FailedSwraData(replyTo, id, new RuntimeException(s"Workflow $id transitioned to state $WorkflowFailed")))
    case Event(RequestComplete((StatusCodes.OK, jsObject: JsObject)), RunningSwraData(replyTo, id)) if jsObject.state == WorkflowAborted =>
      log.info(s"$Tag workflow finished with status '$WorkflowAborted'.")
      requestMetadataOrIssueReply(AbortedSwraData(replyTo, id))
  }

  when (RequestingOutputs) {
    case Event(RequestComplete((StatusCodes.OK, outputs: JsObject)), data: TerminalSwraData) =>
      outputOutputs(outputs)
      requestMetadataOrIssueReply(data)
  }

  when (RequestingMetadata) {
    case Event(RequestComplete((StatusCodes.OK, metadata: JsObject)), data: TerminalSwraData) =>
      outputMetadata(metadata)
      issueReply(data)
  }

  onTransition {
    case NotStarted -> RunningWorkflow => schedulePollRequest()
  }

  whenUnhandled {
    // Handle failures for all failure responses generically.
    case Event(r: WorkflowStoreActor.WorkflowAbortFailed, data) => failAndFinish(r.reason, data)
    case Event(Failure(e), data) => failAndFinish(e, data)
    case Event(Status.Failure(e), data) => failAndFinish(e, data)
    case Event(RequestComplete((_, snap)), data) => failAndFinish(new RuntimeException(s"Unexpected API completion message: $snap"), data)
    case Event((CurrentState(_, _) | Transition(_, _, _)), _) =>
      // ignore uninteresting current state and transition messages
      stay()
    case Event(m, d) =>
      log.warning(s"$Tag: received unexpected message: $m in state ${d.getClass.getSimpleName}")
      stay()
  }

  private def requestMetadataOrIssueReply(newData: TerminalSwraData) = if (metadataOutputPath.isDefined) requestMetadata(newData) else issueReply(newData)
  
  private def requestMetadata(newData: TerminalSwraData): State = {
    val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor), s"MetadataRequest-Workflow-${newData.id}")
    metadataBuilder ! GetSingleWorkflowMetadataAction(newData.id, None, None, expandSubWorkflows = true)
    goto (RequestingMetadata) using newData
  }

  private def schedulePollRequest(): Unit = {
    // -Ywarn-value-discard should stash Cancellable to cancel
    context.system.scheduler.scheduleOnce(backoff.backoffMillis.millis, self, IssuePollRequest)
    ()
  }

  private def requestStatus(id: WorkflowId): Unit = {
    // This requests status via the metadata service rather than instituting an FSM watch on the underlying workflow actor.
    // Cromwell's eventual consistency means it isn't safe to use an FSM transition to a terminal state as the signal for
    // when outputs or metadata have stabilized.
    val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor), s"StatusRequest-Workflow-$id-request-${UUID.randomUUID()}")
    metadataBuilder ! GetStatus(id)
  }

  private def issueSuccessReply(replyTo: ActorRef): State = {
    replyTo.tell(msg = (), sender = self) // Because replyTo ! () is the parameterless call replyTo.!()
    context.stop(self)
    stay()
  }

  private def issueFailureReply(replyTo: ActorRef, e: Throwable): State = {
    replyTo ! Status.Failure(e)
    context.stop(self)
    stay()
  }

  private def issueReply(data: TerminalSwraData) = {
    data match {
      case s: SucceededSwraData => issueSuccessReply(s.replyTo)
      case f: FailedSwraData => issueFailureReply(f.replyTo, f.failure)
      case a: AbortedSwraData => issueSuccessReply(a.replyTo)

    }
  }

  private def failAndFinish(e: Throwable, data: SwraData): State = {
    log.error(e, s"$Tag received Failure message: ${e.getMessage}")
    data match {
      case EmptySwraData =>
        log.error(e, "Cannot issue response. Need a 'replyTo' address to issue the exception response")
        context.stop(self)
        stay()
      case SubmittedSwraData(replyTo) =>
        issueFailureReply(replyTo, e)
      case RunningSwraData(replyTo, _) =>
        issueFailureReply(replyTo, e)
      case c: TerminalSwraData =>
        issueFailureReply(c.replyTo, e)
    }
  }

  /**
    * Outputs the outputs to stdout, and then requests the metadata.
    */
  private def outputOutputs(outputs: JsObject): Unit = {
    println(outputs.prettyPrint)
  }

  private def outputMetadata(metadata: JsObject): Try[Unit] = {
    Try {
      val path = File(metadataOutputPath.get)
      if (path.isDirectory) {
        log.error("Specified metadata path is a directory, should be a file: " + path)
      } else {
        log.info(s"$Tag writing metadata to $path")
        path.createIfNotExists(asDirectory = false, createParents = true).write(metadata.prettyPrint)
      }
    } void
  }
}

object SingleWorkflowRunnerActor {
  def props(source: WorkflowSourceFilesCollection, metadataOutputFile: Option[Path]): Props = {
    Props(new SingleWorkflowRunnerActor(source, metadataOutputFile)).withDispatcher(EngineDispatcher)
  }

  sealed trait RunnerMessage
  // The message to actually run the workflow is made explicit so the non-actor Main can `ask` this actor to do the
  // running and collect a result.
  case object RunWorkflow extends RunnerMessage
  private case object IssuePollRequest extends RunnerMessage

  sealed trait RunnerState
  case object NotStarted extends RunnerState
  case object SubmittedWorkflow extends RunnerState
  case object RunningWorkflow extends RunnerState
  case object RequestingOutputs extends RunnerState
  case object RequestingMetadata extends RunnerState

  sealed trait SwraData
  case object EmptySwraData extends SwraData
  final case class SubmittedSwraData(replyTo: ActorRef) extends SwraData
  final case class RunningSwraData(replyTo: ActorRef, id: WorkflowId) extends SwraData

  sealed trait TerminalSwraData extends SwraData { def replyTo: ActorRef; def terminalState: WorkflowState; def id: WorkflowId }
  final case class SucceededSwraData(replyTo: ActorRef,
                                     id: WorkflowId) extends TerminalSwraData { override val terminalState = WorkflowSucceeded }

  final case class FailedSwraData(replyTo: ActorRef,
                                  id: WorkflowId,
                                  failure: Throwable) extends TerminalSwraData { override val terminalState = WorkflowFailed }

  final case class AbortedSwraData(replyTo: ActorRef,
                                   id: WorkflowId) extends TerminalSwraData { override val terminalState = WorkflowAborted }

  implicit class EnhancedJsObject(val jsObject: JsObject) extends AnyVal {
    def state: WorkflowState = WorkflowState.fromString(jsObject.fields("status").asInstanceOf[JsString].value)
  }

  private val Tag = "SingleWorkflowRunnerActor"
}
