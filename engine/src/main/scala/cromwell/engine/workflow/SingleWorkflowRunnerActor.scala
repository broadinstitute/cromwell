package cromwell.engine.workflow

import java.util.UUID

import akka.actor.FSM.{CurrentState, Transition}
import akka.actor._
import akka.stream.ActorMaterializer
import cats.instances.try_._
import cats.syntax.functor._
import com.typesafe.config.Config
import common.util.VersionUtil
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.core.abort.WorkflowAbortFailureResponse
import cromwell.core.actor.BatchActor.QueueWeight
import cromwell.core.path.Path
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.engine.CromwellTerminator
import cromwell.engine.workflow.SingleWorkflowRunnerActor._
import cromwell.engine.workflow.WorkflowManagerActor.{PreventNewWorkflowsFromStarting, RetrieveNewWorkflows}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.SubmitWorkflow
import cromwell.engine.workflow.workflowstore.{InMemoryWorkflowStore, WorkflowStoreSubmitActor}
import cromwell.jobstore.EmptyJobStoreActor
import cromwell.server.CromwellRootActor
import cromwell.services.metadata.MetadataService.{GetSingleWorkflowMetadataAction, GetStatus, ListenToMetadataWriteActor, WorkflowOutputs}
import cromwell.subworkflowstore.EmptySubWorkflowStoreActor
import cromwell.webservice.metadata.MetadataBuilderActor
import cromwell.webservice.metadata.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse}
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Try}

/**
 * Designed explicitly for the use case of the 'run' functionality in Main. This Actor will start a workflow,
 * print out the outputs when complete and reply with a result.
 */
class SingleWorkflowRunnerActor(source: WorkflowSourceFilesCollection,
                                metadataOutputPath: Option[Path],
                                terminator: CromwellTerminator,
                                gracefulShutdown: Boolean,
                                abortJobsOnTerminate: Boolean,
                                config: Config
                                )(implicit materializer: ActorMaterializer)
  extends CromwellRootActor(terminator, gracefulShutdown, abortJobsOnTerminate, false, config)
    with LoggingFSM[RunnerState, SwraData] {

  import SingleWorkflowRunnerActor._
  private val backoff = SimpleExponentialBackoff(1 second, 1 minute, 1.2)

  override lazy val workflowStore = new InMemoryWorkflowStore()
  override lazy val jobStoreActor = context.actorOf(EmptyJobStoreActor.props, "JobStoreActor")
  override lazy val subWorkflowStoreActor = context.actorOf(EmptySubWorkflowStoreActor.props, "SubWorkflowStoreActor")

  log.info("{}: Version {}", Tag, VersionUtil.getVersion("cromwell-engine"))
  startWith(NotStarted, EmptySwraData)

  when (NotStarted) {
    case Event(RunWorkflow, EmptySwraData) =>
      log.info(s"$Tag: Submitting workflow")
      workflowStoreActor ! SubmitWorkflow(source)
      goto(SubmittedWorkflow) using SubmittedSwraData(sender())
  }

  when (SubmittedWorkflow) {
    case Event(WorkflowStoreSubmitActor.WorkflowSubmittedToStore(id, WorkflowSubmitted), SubmittedSwraData(replyTo)) =>
      log.info(s"$Tag: Workflow submitted UUID($id)")
      // Since we only have a single workflow, force the WorkflowManagerActor's hand in case the polling rate is long
      workflowManagerActor ! RetrieveNewWorkflows
      // After that - prevent the WMA from trying to start new workflows
      workflowManagerActor ! PreventNewWorkflowsFromStarting
      schedulePollRequest()
      goto(RunningWorkflow) using RunningSwraData(replyTo, id)
  }

  when (RunningWorkflow) {
    case Event(IssuePollRequest, RunningSwraData(_, id)) =>
      requestStatus(id)
      stay()
    case Event(BuiltMetadataResponse(jsObject: JsObject), RunningSwraData(_, _)) if !jsObject.state.isTerminal =>
      schedulePollRequest()
      stay()
    case Event(BuiltMetadataResponse(jsObject: JsObject), RunningSwraData(replyTo, id)) if jsObject.state == WorkflowSucceeded =>
      log.info(s"$Tag workflow finished with status '$WorkflowSucceeded'.")
      serviceRegistryActor ! ListenToMetadataWriteActor
      goto(WaitingForFlushedMetadata) using SucceededSwraData(replyTo, id)
    case Event(BuiltMetadataResponse(jsObject: JsObject), RunningSwraData(replyTo, id)) if jsObject.state == WorkflowFailed =>
      log.info(s"$Tag workflow finished with status '$WorkflowFailed'.")
      serviceRegistryActor ! ListenToMetadataWriteActor
      goto(WaitingForFlushedMetadata) using FailedSwraData(replyTo, id, new RuntimeException(s"Workflow $id transitioned to state $WorkflowFailed"))
    case Event(BuiltMetadataResponse(jsObject: JsObject), RunningSwraData(replyTo, id)) if jsObject.state == WorkflowAborted =>
      log.info(s"$Tag workflow finished with status '$WorkflowAborted'.")
      serviceRegistryActor ! ListenToMetadataWriteActor
      goto(WaitingForFlushedMetadata) using AbortedSwraData(replyTo, id)
  }
  
  when (WaitingForFlushedMetadata) {
    case Event(QueueWeight(weight), _) if weight > 0 => stay()
    case Event(QueueWeight(_), data: SucceededSwraData) =>
      val metadataBuilder = context.actorOf(MetadataBuilderActor.props(serviceRegistryActor),
        s"CompleteRequest-Workflow-${data.id}-request-${UUID.randomUUID()}")
      metadataBuilder ! WorkflowOutputs(data.id)
      goto(RequestingOutputs)
    case Event(QueueWeight(_), data : TerminalSwraData) =>
      requestMetadataOrIssueReply(data)
  }

  when (RequestingOutputs) {
    case Event(BuiltMetadataResponse(outputs: JsObject), data: TerminalSwraData) =>
      outputOutputs(outputs)
      requestMetadataOrIssueReply(data)
  }

  when (RequestingMetadata) {
    case Event(BuiltMetadataResponse(metadata: JsObject), data: TerminalSwraData) =>
      outputMetadata(metadata)
      issueReply(data)
  }

  onTransition {
    case NotStarted -> RunningWorkflow => schedulePollRequest()
  }

  whenUnhandled {
    // Handle failures for all failure responses generically.
    case Event(r: WorkflowAbortFailureResponse, data) => failAndFinish(r.failure, data)
    case Event(Failure(e), data) => failAndFinish(e, data)
    case Event(Status.Failure(e), data) => failAndFinish(e, data)
    case Event(FailedMetadataResponse(e), data) => failAndFinish(e, data)
    case Event(CurrentState(_, _) | Transition(_, _, _), _) =>
      // ignore uninteresting current state and transition messages
      stay()
    case Event(akka.Done, _) =>
      // This actor sends the WMA a `PreventNewWorkflowsFromStarting` message to which the WMA will respond with
      // an `akka.Done` message, so don't act surprised when that happens.
      stay()
    case Event(m, _) =>
      log.warning(s"$Tag: received unexpected message: $m of type ${m.getClass.getCanonicalName} in state $stateName")
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
    done()
    stay()
  }

  private def issueFailureReply(replyTo: ActorRef, e: Throwable): State = {
    replyTo ! Status.Failure(e)
    done()
    stay()
  }
  
  private [workflow] def done() = {}

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
        done()
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
      val path = metadataOutputPath.get
      if (path.isDirectory) {
        log.error("Specified metadata path is a directory, should be a file: " + path)
      } else {
        log.info(s"$Tag writing metadata to $path")
        path.createIfNotExists(createParents = true).write(metadata.prettyPrint)
      }
    } void
  }
}

object SingleWorkflowRunnerActor {
  def props(source: WorkflowSourceFilesCollection,
            metadataOutputFile: Option[Path],
            terminator: CromwellTerminator,
            gracefulShutdown: Boolean,
            abortJobsOnTerminate: Boolean,
            config: Config)
           (implicit materializer: ActorMaterializer): Props = {
    Props(
      new SingleWorkflowRunnerActor(
        source = source,
        metadataOutputPath = metadataOutputFile,
        terminator = terminator,
        gracefulShutdown = gracefulShutdown,
        abortJobsOnTerminate = abortJobsOnTerminate,
        config = config
      )
    ).withDispatcher(EngineDispatcher)
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
  case object WaitingForFlushedMetadata extends RunnerState
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
    def state: WorkflowState = WorkflowState.withName(jsObject.fields("status").asInstanceOf[JsString].value)
  }

  private val Tag = "SingleWorkflowRunnerActor"
}
