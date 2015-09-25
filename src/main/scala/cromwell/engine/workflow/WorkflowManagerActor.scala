package cromwell.engine.workflow

import akka.actor.FSM.{SubscribeTransitionCallBack, Transition}
import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.pipe
import com.typesafe.config.ConfigFactory
import cromwell.binding
import cromwell.binding._
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.{Backend, CallMetadata, StdoutStderr}
import cromwell.engine.db.slick._
import cromwell.engine.db.{DataAccess, ExecutionDatabaseKey}
import cromwell.engine.workflow.WorkflowActor.{Restart, Start}
import cromwell.util.WriteOnceStore
import cromwell.webservice.WorkflowMetadataResponse
import org.joda.time.DateTime
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.io.Source
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WorkflowManagerActor {
  class WorkflowNotFoundException(message: String) extends RuntimeException(message)
  class CallNotFoundException(message: String) extends RuntimeException(message)

  sealed trait WorkflowManagerActorMessage
  case class SubmitWorkflow(source: WorkflowSourceFiles) extends WorkflowManagerActorMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerActorMessage
  case class CallOutputs(id: WorkflowId, callFqn: FullyQualifiedName) extends WorkflowManagerActorMessage
  case class CallStdoutStderr(id: WorkflowId, callFqn: FullyQualifiedName) extends WorkflowManagerActorMessage
  case class WorkflowStdoutStderr(id: WorkflowId) extends WorkflowManagerActorMessage
  case object Shutdown extends WorkflowManagerActorMessage
  case class SubscribeToWorkflow(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowAbort(id: WorkflowId) extends WorkflowManagerActorMessage
  final case class WorkflowMetadata(id: WorkflowId) extends WorkflowManagerActorMessage

  def props(backend: Backend): Props = Props(new WorkflowManagerActor(backend))

  lazy val BackendInstance = Backend.from(ConfigFactory.load.getConfig("backend"))
  lazy val BackendType = BackendInstance.backendType
}

/**
 * Responses to messages:
 * SubmitWorkflow: Returns a `Future[WorkflowId]`
 * WorkflowStatus: Returns a `Future[Option[WorkflowState]]`
 * WorkflowOutputs: Returns a `Future[Option[binding.WorkflowOutputs]]` aka `Future[Option[Map[String, WdlValue]]`
 *
 */
class WorkflowManagerActor(backend: Backend) extends Actor with CromwellActor {
  import WorkflowManagerActor._
  private val log = Logging(context.system, this)
  private val tag = "WorkflowManagerActor"

  type WorkflowActorRef = ActorRef

  private val workflowStore = new WriteOnceStore[WorkflowId, WorkflowActorRef]

  override def preStart() {
    restartIncompleteWorkflows()
  }

  def receive = LoggingReceive {
    case SubmitWorkflow(source) =>
      submitWorkflow(source, maybeWorkflowId = None) pipeTo sender
    case WorkflowStatus(id) => DataAccess.instance.getWorkflowState(id) pipeTo sender
    case WorkflowAbort(id) =>
      workflowStore.toMap.get(id) match {
        case Some(x) =>
          x ! WorkflowActor.AbortWorkflow
          sender ! Some(WorkflowAborting)
        case None => sender ! None
      }
    case Shutdown => context.system.shutdown()
    case WorkflowOutputs(id) => workflowOutputs(id) pipeTo sender
    case CallOutputs(workflowId, callName) => callOutputs(workflowId, callName) pipeTo sender
    case CallStdoutStderr(workflowId, callName) => callStdoutStderr(workflowId, callName) pipeTo sender
    case WorkflowStdoutStderr(workflowId) => workflowStdoutStderr(workflowId) pipeTo sender
    case WorkflowMetadata(workflowId) => workflowMetadata(workflowId) pipeTo sender
    case Transition(actor, oldState, newState: WorkflowState) => updateWorkflowState(actor, newState)
    case SubscribeToWorkflow(id) =>
      //  NOTE: This fails silently. Currently we're ok w/ this, but you might not be in the future
      workflowStore.toMap.get(id) foreach {_ ! SubscribeTransitionCallBack(sender())}
  }

  /**
   * Returns a `Future[Any]` which will be failed if there is no workflow with the specified id.
   */
  private def assertWorkflowExistence(id: WorkflowId): Future[Any] = {
    // Confirm the workflow exists by querying its state.  If no state is found the workflow doesn't exist.
    DataAccess.instance.getWorkflowState(id) map {
      case None => throw new WorkflowNotFoundException(s"Workflow '$id' not found")
      case _ =>
    }
  }

  private def assertCallExistence(id: WorkflowId, callFqn: FullyQualifiedName): Future[Any] = {
    DataAccess.instance.getExecutionStatus(id, ExecutionDatabaseKey(callFqn, None)) map {
      case None => throw new CallNotFoundException(s"Call '$callFqn' not found in workflow '$id'.")
      case _ =>
    }
  }

  /**
   * Retrieve the entries that produce stdout and stderr.
   */
  private def getCallLogKeys(id: WorkflowId, callFqn: FullyQualifiedName): Future[Seq[ExecutionDatabaseKey]] = {
    DataAccess.instance.getExecutionStatuses(id, callFqn) map {
      case map if map.isEmpty => throw new CallNotFoundException(s"Call '$callFqn' not found in workflow '$id'.")
      case entries =>
        val callKeys = entries.keys filterNot isCollector(entries.keys)
        callKeys.toSeq.sortBy(_.index)
    }
  }

  private def workflowOutputs(id: WorkflowId): Future[binding.WorkflowOutputs] = {
    for {
      _ <- assertWorkflowExistence(id)
      outputs <- DataAccess.instance.getOutputs(id)
    } yield {
      SymbolStoreEntry.toWorkflowOutputs(outputs)
    }
  }

  private def callOutputs(workflowId: WorkflowId, callFqn: String): Future[binding.CallOutputs] = {
    for {
      _ <- assertWorkflowExistence(workflowId)
      _ <- assertCallExistence(workflowId, callFqn)
      outputs <- DataAccess.instance.getOutputs(workflowId, ExecutionDatabaseKey(callFqn, None))
    } yield {
      SymbolStoreEntry.toCallOutputs(outputs)
    }
  }

  private def assertCallFqnWellFormed(descriptor: WorkflowDescriptor, callFqn: FullyQualifiedName): Try[String] = {
    descriptor.namespace.resolve(callFqn) match {
      case Some(c: Call) => Success(c.name)
      case _ => Failure(new UnsupportedOperationException("Expected a fully qualified name to have at least two parts"))
    }
  }

  private def isCollector(entries: Iterable[ExecutionDatabaseKey])(key: ExecutionDatabaseKey) = {
    key.index.isEmpty &&
      (entries exists { e =>
        (e.fqn == key.fqn) && e.index.isDefined
      })
  }

  private def hasLogs(entries: Iterable[ExecutionDatabaseKey])(key: ExecutionDatabaseKey) = {
    !key.fqn.isScatter && !isCollector(entries)(key)
  }

  private def callStdoutStderr(workflowId: WorkflowId, callFqn: String): Future[Any] = {
    for {
        _ <- assertWorkflowExistence(workflowId)
        descriptor <- DataAccess.instance.getWorkflow(workflowId)
        callName <- Future.fromTry(assertCallFqnWellFormed(descriptor, callFqn))
        callLogKeys <- getCallLogKeys(workflowId, callFqn)
        callStandardOutput <- Future.successful(callLogKeys map { key => backend.stdoutStderr(descriptor, callName, key.index) })
      } yield callStandardOutput
  }

  private def workflowStdoutStderr(workflowId: WorkflowId): Future[Map[FullyQualifiedName, Seq[StdoutStderr]]] = {
    def logMapFromStatusMap(descriptor: WorkflowDescriptor, statusMap: Map[ExecutionDatabaseKey, ExecutionStatus]): Try[Map[FullyQualifiedName, Seq[StdoutStderr]]] = {
      Try {
        val sortedMap = statusMap.toSeq.sortBy(_._1.index)
        val callsToPaths = for {
          (key, status) <- sortedMap if status.isTerminal && hasLogs(statusMap.keys)(key)
          callName = assertCallFqnWellFormed(descriptor, key.fqn).get
          callStandardOutput = backend.stdoutStderr(descriptor, callName, key.index)
        } yield key.fqn -> callStandardOutput

        callsToPaths groupBy { _._1 } mapValues { v => v map { _._2 } }
      }
    }

    for {
      _ <- assertWorkflowExistence(workflowId)
      descriptor <- DataAccess.instance.getWorkflow(workflowId)
      callToStatusMap <- DataAccess.instance.getExecutionStatuses(workflowId)
      x = callToStatusMap mapValues { _.executionStatus }
      callToLogsMap <- Future.fromTry(logMapFromStatusMap(descriptor, callToStatusMap mapValues { _.executionStatus }))
    } yield callToLogsMap
  }

  private def buildCallMetadata(executions: Traversable[Execution],
                                standardStreamsMap: Map[FullyQualifiedName, Seq[StdoutStderr]],
                                callInputs: Traversable[Symbol],
                                callOutputs: Traversable[Symbol],
                                jobInfos: Map[ExecutionDatabaseKey, Any]): Map[FullyQualifiedName, Seq[CallMetadata]] = {

    // add straightforward queries for call inputs and outputs
    // navigate the keys of the map and indices of their arrays to effectively unify the data
    // into CallMetadata.
    // figure out how to get job IDs and maybe on which backend the call executed
    standardStreamsMap map { case (key, seqOfStreams) =>
      key -> seqOfStreams.map { streams =>
        CallMetadata(
          inputs = Map("input_key" -> "input_value"),
          status = "UnknownStatus",
          backend = Option("UnknownBackend"),
          outputs = Option(Map("output_key" -> "output_value")),
          start = Option(new DateTime()),
          end = Option(new DateTime()),
          jobId = Option("COMPLETELY-MADE-UP-ID"),
          rc = Option(0),
          stdout = Option(streams.stdout),
          stderr = Option(streams.stderr))
      }
    }
  }

  private def buildWorkflowMetadata(workflowExecution: WorkflowExecution,
                                    workflowExecutionAux: WorkflowExecutionAux,
                                    workflowOutputs: binding.WorkflowOutputs,
                                    callMetadata: Map[FullyQualifiedName, Seq[CallMetadata]]): WorkflowMetadataResponse = {

    val startDate = new DateTime(workflowExecution.startDt)
    val endDate = workflowExecution.endDt map { new DateTime(_) }
    val outputs = Option(workflowOutputs mapValues { _.valueString })
    val workflowInputsString = Source.fromInputStream(workflowExecutionAux.jsonInputs.getAsciiStream).mkString
    // The casting here looks rough but should be safe if the workflow has gotten through submission?
    val workflowInputsMap = workflowInputsString.parseJson.asInstanceOf[JsObject].fields map {
      case (k, v) => k -> v.asInstanceOf[JsString].value
    }

    WorkflowMetadataResponse(
      id = workflowExecution.workflowExecutionUuid.toString,
      status = workflowExecution.status,
      // We currently do not make a distinction between the submission and start dates of a workflow, but it's
      // possible at least theoretically that a workflow might not begin to execute immediately upon submission.
      submission = startDate,
      start = Option(startDate),
      end = endDate,
      inputs = workflowInputsMap,
      outputs = outputs,
      calls = callMetadata)
  }

  private def workflowMetadata(id: WorkflowId): Future[WorkflowMetadataResponse] = {
    for {
      workflowExecution <- DataAccess.instance.getWorkflowExecution(id)
      workflowOutputs <- workflowOutputs(id)
      // The workflow has been persisted in the DB so we know the workflowExecutionId must be non-null,
      // so the .get on the Option is safe.
      workflowExecutionAux <- DataAccess.instance.getWorkflowExecutionAux(workflowExecution.workflowExecutionId.get)
      callStandardStreamsMap <- workflowStdoutStderr(id)
      executions <- DataAccess.instance.getExecutions(workflowExecution.workflowExecutionId.get)
      callInputs <- DataAccess.instance.getAllInputs(id)
      callOutputs <- DataAccess.instance.getAllOutputs(id)
      jesJobs <- DataAccess.instance.jesJobInfo(id)
      localJobs <- DataAccess.instance.localJobInfo(id)
      sgeJobs <- DataAccess.instance.sgeJobInfo(id)

      callMetadata = buildCallMetadata(executions, callStandardStreamsMap, callInputs, callOutputs, jesJobs ++ localJobs ++ sgeJobs)
      workflowMetadata = buildWorkflowMetadata(workflowExecution, workflowExecutionAux, workflowOutputs, callMetadata)

    } yield workflowMetadata
  }

  private def submitWorkflow(source: WorkflowSourceFiles,
                             maybeWorkflowId: Option[WorkflowId]): Future[WorkflowId] = {
    val workflowId: WorkflowId = maybeWorkflowId.getOrElse(WorkflowId.randomId())
    log.info(s"$tag submitWorkflow input id = $maybeWorkflowId, effective id = $workflowId")
    val futureId = for {
      descriptor <- Future.fromTry(Try(new WorkflowDescriptor(workflowId, source)))
      workflowActor = context.actorOf(WorkflowActor.props(descriptor, backend), s"WorkflowActor-$workflowId")
      _ <- Future.fromTry(workflowStore.insert(workflowId, workflowActor))
    } yield {
      val isRestart = maybeWorkflowId.isDefined
      workflowActor ! (if (isRestart) Restart else Start)
      workflowActor ! SubscribeTransitionCallBack(self)
      workflowId
    }

    futureId onFailure {
      case e =>
        val messageOrBlank = Option(e.getMessage).mkString
        log.error(e, s"$tag: Workflow failed submission: " + messageOrBlank)
    }

    futureId
  }

  private def restartWorkflow(restartableWorkflow: RestartableWorkflow): Unit =
    submitWorkflow(restartableWorkflow.source, Option(restartableWorkflow.id))

  private def updateWorkflowState(workflow: WorkflowActorRef, state: WorkflowState): Future[Unit] = {
    val id = idByWorkflow(workflow)
    DataAccess.instance.updateWorkflowState(id, state)
  }

  private def idByWorkflow(workflow: WorkflowActorRef): WorkflowId = {
    workflowStore.toMap collectFirst { case (k, v) if v == workflow => k } get
  }

  private def restartIncompleteWorkflows(): Unit = {
    type QuantifierAndPlural = (String, String)
    def pluralize(num: Int): QuantifierAndPlural = {
      (if (num == 0) "no" else num.toString, if (num == 1) "" else "s")
    }

    def buildRestartableWorkflows(workflowDescriptors: Traversable[WorkflowDescriptor]): Seq[RestartableWorkflow] = {
      (for {
        workflowDescriptor <- workflowDescriptors
      } yield RestartableWorkflow(workflowDescriptor.id, workflowDescriptor.sourceFiles)).toSeq
    }

    val result = for {
      workflowDescriptors <- DataAccess.instance.getWorkflowsByState(Seq(WorkflowSubmitted, WorkflowRunning))
      restartableWorkflows = buildRestartableWorkflows(workflowDescriptors)
      _ <- backend.handleCallRestarts(restartableWorkflows)
    } yield {
        val num = restartableWorkflows.length
        val (displayNum, plural) = pluralize(num)
        log.info(s"$tag Found $displayNum workflow$plural to restart.")

        if (num > 0) {
          val ids = restartableWorkflows.map { _.id.toString }.sorted
          log.info(s"$tag Restarting workflow ID$plural: " + ids.mkString(", "))
        }

        restartableWorkflows foreach restartWorkflow
    }

    result recover {
      case e: Throwable => log.error(e, e.getMessage)
    }
  }
}
