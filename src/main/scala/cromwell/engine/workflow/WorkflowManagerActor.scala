package cromwell.engine.workflow

import java.util.UUID

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor.{Actor, ActorRef, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.pipe
import com.typesafe.config.ConfigFactory
import cromwell.binding
import cromwell.binding._
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.{Backend, StdoutStderr}
import cromwell.engine.db.{ExecutionDatabaseKey, DataAccess}
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.workflow.WorkflowActor.{Restart, Start}
import cromwell.util.WriteOnceStore
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}


object WorkflowManagerActor {
  class WorkflowNotFoundException(message: String) extends RuntimeException(message)
  class CallNotFoundException(message: String) extends RuntimeException(message)

  sealed trait WorkflowManagerActorMessage
  case class SubmitWorkflow(wdlSource: WdlSource, wdlJson: WdlJson, inputs: binding.WorkflowRawInputs) extends WorkflowManagerActorMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerActorMessage
  case class CallOutputs(id: WorkflowId, callFqn: FullyQualifiedName) extends WorkflowManagerActorMessage
  case class CallStdoutStderr(id: WorkflowId, callFqn: FullyQualifiedName) extends WorkflowManagerActorMessage
  case class WorkflowStdoutStderr(id: WorkflowId) extends WorkflowManagerActorMessage
  case object Shutdown extends WorkflowManagerActorMessage
  case class SubscribeToWorkflow(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowAbort(id: WorkflowId) extends WorkflowManagerActorMessage

  def props(dataAccess: DataAccess, backend: Backend): Props = Props(new WorkflowManagerActor(dataAccess, backend))

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
class WorkflowManagerActor(dataAccess: DataAccess, backend: Backend) extends Actor with CromwellActor {
  import WorkflowManagerActor._
  private val log = Logging(context.system, this)
  private val tag = "WorkflowManagerActor"

  type WorkflowActorRef = ActorRef

  private val workflowStore = new WriteOnceStore[WorkflowId, WorkflowActorRef]

  override def preStart() {
    restartIncompleteWorkflows()
  }

  def receive = LoggingReceive {
    case SubmitWorkflow(wdlSource, wdlJson, inputs) =>
      submitWorkflow(wdlSource, wdlJson, inputs, maybeWorkflowId = None) pipeTo sender
    case WorkflowStatus(id) => dataAccess.getWorkflowState(id) pipeTo sender
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
    case CurrentState(actor, state: WorkflowState) => updateWorkflowState(actor, state)
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
    dataAccess.getWorkflowState(id) map {
      case None => throw new WorkflowNotFoundException(s"Workflow '$id' not found")
      case _ =>
    }
  }

  private def assertCallExistence(id: WorkflowId, callFqn: FullyQualifiedName): Future[Any] = {
    dataAccess.getExecutionStatus(id, ExecutionDatabaseKey(callFqn, None)) map {
      case None => throw new CallNotFoundException(s"Call '$callFqn' not found in workflow '$id'.")
      case _ =>
    }
  }

  private def workflowOutputs(id: WorkflowId): Future[binding.WorkflowOutputs] = {
    for {
      _ <- assertWorkflowExistence(id)
      outputs <- dataAccess.getOutputs(id)
    } yield {
      SymbolStoreEntry.toWorkflowOutputs(outputs)
    }
  }

  private def callOutputs(workflowId: WorkflowId, callFqn: String): Future[binding.CallOutputs] = {
    for {
      _ <- assertWorkflowExistence(workflowId)
      _ <- assertCallExistence(workflowId, callFqn)
      outputs <- dataAccess.getOutputs(workflowId, ExecutionDatabaseKey(callFqn, None))
    } yield {
      SymbolStoreEntry.toCallOutputs(outputs)
    }
  }

  /* Return value here is a tuple: (workflow name, call name)
   * TODO: This should use a more sophisticated way of resolving FQNs.  See DSDEEPB-986
   */
  private def assertCallFqnWellFormed(callFqn: FullyQualifiedName): Try[(String, String)] = {
    callFqn.split("\\.").toSeq match {
      case s: Seq[String] if s.size >= 2 => Success(s.head, s.last)
      case _ => Failure(new UnsupportedOperationException("Expected a fully qualified name to have at least two parts"))
    }
  }

  private def callStdoutStderr(workflowId: WorkflowId, callFqn: String): Future[StdoutStderr] = {
    for {
      _ <- assertWorkflowExistence(workflowId)
      _ <- assertCallExistence(workflowId, callFqn)
      (wf, call) <- Future.fromTry(assertCallFqnWellFormed(callFqn))
      callStandardOutput <- Future.successful(backend.stdoutStderr(workflowId, wf, call))
    } yield callStandardOutput
  }

  private def workflowStdoutStderr(workflowId: WorkflowId): Future[Map[FullyQualifiedName, StdoutStderr]] = {
    def logMapFromStatusMap(statusMap: Map[FullyQualifiedName, ExecutionStatus]): Try[Map[FullyQualifiedName, StdoutStderr]] = {
      Try {
        val callsToPaths = for {
          (call, status) <- statusMap.toSeq
          if Set(ExecutionStatus.Done, ExecutionStatus.Failed).contains(status)
          (wf, callLqn) = assertCallFqnWellFormed(call).get
          callStandardOutput = backend.stdoutStderr(workflowId, wf, callLqn)
        } yield call -> callStandardOutput
        callsToPaths.toMap
      }
    }

    for {
      _ <- assertWorkflowExistence(workflowId)
      callToStatusMap <- dataAccess.getExecutionStatuses(workflowId)
      callToLogsMap <- Future.fromTry(logMapFromStatusMap(callToStatusMap map { case (k, v) => (k.fqn, v) }))
    } yield callToLogsMap
  }

  private def submitWorkflow(wdlSource: WdlSource, wdlJson: WdlJson, inputs: WorkflowRawInputs,
                             maybeWorkflowId: Option[WorkflowId]): Future[WorkflowId] = {
    val workflowId: WorkflowId = maybeWorkflowId.getOrElse(WorkflowId.randomId())
    log.info(s"$tag submitWorkflow input id = $maybeWorkflowId, effective id = $workflowId")
    val futureId = for {
      eventualNamespace <- Future(NamespaceWithWorkflow.load(wdlSource, BackendType))
      coercedInputs <- Future.fromTry(eventualNamespace.coerceRawInputs(inputs))
      declarations <- Future.fromTry(eventualNamespace.staticDeclarationsRecursive(coercedInputs))
      inputs = coercedInputs ++ declarations
      descriptor = new WorkflowDescriptor(workflowId, eventualNamespace, wdlSource, wdlJson, inputs)
      workflowActor = context.actorOf(WorkflowActor.props(descriptor, backend, dataAccess), s"WorkflowActor-$workflowId")
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

  private def restartWorkflow(restartableWorkflow: RestartableWorkflow): Unit = {
    submitWorkflow(restartableWorkflow.source, restartableWorkflow.json, restartableWorkflow.inputs, Option(restartableWorkflow.id))
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
          log.error(s"$tag Error restarting workflow ${workflowInfo.workflowId}: expected JSON inputs, got '$x'.")
          None
      }
    }

    type QuantifierAndPlural = (String, String)
    def pluralize(num: Int): QuantifierAndPlural = {
      (if (num == 0) "no" else num.toString, if (num == 1) "" else "s")
    }

    def buildRestartableWorkflows(workflowInfos: Traversable[WorkflowInfo]): Seq[RestartableWorkflow] = {
      (for {
        workflowInfo <- workflowInfos
        jsonInputs = clobToJsonInputs(workflowInfo)
        if jsonInputs.isDefined
      } yield RestartableWorkflow(workflowInfo.workflowId, workflowInfo.wdlSource, workflowInfo.wdlJson, jsonInputs.get)).toSeq
    }

    val result = for {
      workflowInfos <- dataAccess.getWorkflowsByState(Seq(WorkflowSubmitted, WorkflowRunning))
      restartableWorkflows = buildRestartableWorkflows(workflowInfos)
      _ <- backend.handleCallRestarts(restartableWorkflows, dataAccess)
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
