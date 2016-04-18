package cromwell.engine.workflow

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor._
import akka.event.Logging
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.{ExecutionDatabaseKey, ExecutionInfosByExecution}
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorFailure, MaterializeWorkflowDescriptorSuccess}
import cromwell.engine.workflow.ShadowWorkflowActor._
import cromwell.engine.workflow.WorkflowActor.{Restart, Start}
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.util.PromiseActor
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice._
import cromwell.{core, engine}
import lenthall.config.ScalaConfig.EnhancedScalaConfig
import wdl4s._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Future, Promise}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WorkflowManagerActor {
  class WorkflowNotFoundException(message: String) extends RuntimeException(message)
  class CallNotFoundException(message: String) extends RuntimeException(message)

  type WorkflowActorRef = ActorRef
  type WorkflowIdToActorRef = (WorkflowId, WorkflowActorRef)

  sealed trait WorkflowManagerActorMessage

  case class SubmitWorkflow(source: WorkflowSourceFiles) extends WorkflowManagerActorMessage
  case class WorkflowActorSubmitSuccess(replyTo: Option[ActorRef], id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowActorSubmitFailure(replyTo: Option[ActorRef], failure: Throwable) extends WorkflowManagerActorMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowQuery(parameters: Seq[(String, String)]) extends WorkflowManagerActorMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerActorMessage
  case class CallOutputs(id: WorkflowId, callFqn: FullyQualifiedName) extends WorkflowManagerActorMessage
  case class CallStdoutStderr(id: WorkflowId, callFqn: FullyQualifiedName) extends WorkflowManagerActorMessage
  case class WorkflowStdoutStderr(id: WorkflowId) extends WorkflowManagerActorMessage
  case class SubscribeToWorkflow(id: WorkflowId) extends WorkflowManagerActorMessage
  case class WorkflowAbort(id: WorkflowId) extends WorkflowManagerActorMessage
  final case class WorkflowMetadata(id: WorkflowId) extends WorkflowManagerActorMessage
  final case class RestartWorkflows(workflows: Seq[WorkflowDescriptor]) extends WorkflowManagerActorMessage
  final case class CallCaching(id: WorkflowId, parameters: QueryParameters, call: Option[String]) extends WorkflowManagerActorMessage
  private final case class AddEntryToWorkflowManagerData(entry: WorkflowIdToActorRef) extends WorkflowManagerActorMessage
  case object AbortAllWorkflows extends WorkflowManagerActorMessage

  def props(shadowMode: Boolean): Props =
    if (shadowMode) Props(new ShadowWorkflowManagerActor())
    else Props(new WorkflowManagerActor())

  // FIXME hack to deal with one class of "singularity" where Cromwell isn't smart enough to launch only
  // as much work as can reasonably be handled.
  // How long to delay between restarting each workflow that needs to be restarted.  Attempting to
  // restart 500 workflows at exactly the same time crushes the database connection pool.
  lazy val RestartDelay = 200 milliseconds

  sealed trait WorkflowManagerState
  case object Running extends WorkflowManagerState
  case object Aborting extends WorkflowManagerState
  case object Done extends WorkflowManagerState

  case class WorkflowManagerData(workflows: Map[WorkflowId, WorkflowActorRef]) {
    def add(entry: (WorkflowId, WorkflowActorRef)): WorkflowManagerData = this.copy(workflows = workflows + entry)
    def remove(id: WorkflowId): WorkflowManagerData = this.copy(workflows = workflows - id)
    def remove(actor: WorkflowActorRef): WorkflowManagerData = {
      val workflowId = workflows.collectFirst { case (id, a) if a == actor => id }
      // If the ID was found in the lookup return a modified copy of the state data, otherwise just return
      // the same state data.
      workflowId map remove getOrElse this
    }
  }

  lazy val defaultConfig = ConfigFactory.load
}

/**
 * Responses to messages:
 * SubmitWorkflow: Returns a `Future[WorkflowId]`
 * WorkflowStatus: Returns a `Future[Option[WorkflowState]]`
 * WorkflowOutputs: Returns a `Future[Option[binding.WorkflowOutputs]]` aka `Future[Option[Map[String, WdlValue]]`
 *
 */
class WorkflowManagerActor(config: Config)
  extends LoggingFSM[WorkflowManagerState, WorkflowManagerData] with CromwellActor {

  def this() = this(WorkflowManagerActor.defaultConfig)
  implicit val actorSystem = context.system

  private val logger = Logging(context.system, this)
  private val tag = "WorkflowManagerActor"

  private val donePromise = Promise[Unit]()

  override def preStart() {
    addShutdownHook()
    restartIncompleteWorkflows()
  }

  private def addShutdownHook(): Unit = {
    // Only abort jobs on SIGINT if the config explicitly sets backend.abortJobsOnTerminate = true.
    val abortJobsOnTerminate =
      config.getConfig("backend").getBooleanOr("abortJobsOnTerminate", default = false)

    if (abortJobsOnTerminate) {
      sys.addShutdownHook {
        logger.info(s"$tag: Received shutdown signal. Aborting all running workflows...")
        self ! AbortAllWorkflows
        Await.ready(donePromise.future, Duration.Inf)
      }
    }
  }

  private def reply[A](id: WorkflowId, fa: Future[A],
                       onSuccess: (WorkflowId, A) => WorkflowManagerSuccessResponse,
                       onFailure: (WorkflowId, Throwable) => WorkflowManagerFailureResponse): Unit = {
    val sndr = sender()
    fa onComplete {
      case Success(a) => sndr ! onSuccess(id, a)
      case Failure(e) => sndr ! onFailure(id, e)
    }
  }

  private def reply[A, B](id: WorkflowId, a: A, fb: Future[B],
                          onSuccess: (WorkflowId, A, B) => WorkflowManagerSuccessResponse,
                          onFailure: (WorkflowId, A, Throwable) => WorkflowManagerFailureResponse): Unit = {
    val sndr = sender()
    fb onComplete {
      case Success(b) => sndr ! onSuccess(id, a, b)
      case Failure(e) => sndr ! onFailure(id, a, e)
    }
  }

  private def replyToAbort(id: WorkflowId): Unit = {
    // Access the mutable state here synchronously, not in the callback.
    val workflowActor = stateData.workflows.get(id)
    val sndr = sender()
    globalDataAccess.getWorkflowState(id) onComplete {
      case Success(Some(s)) if s.isTerminal =>
        sndr ! WorkflowManagerAbortFailure(id, new IllegalStateException(s"Workflow $id is in terminal state $s and cannot be aborted."))
      case Success(Some(s)) =>
        workflowActor foreach { _ ! WorkflowActor.AbortWorkflow }
        sndr ! WorkflowManagerAbortSuccess(id)
      case Failure(x) => sender ! WorkflowManagerAbortFailure(id, x)
    }
  }

  private def replyToQuery(rawParameters: Seq[(String, String)]): Unit = {
    val sndr = sender()
    query(rawParameters) onComplete {
      case Success(r) => sndr ! WorkflowManagerQuerySuccess(r)
      case Failure(e) => sndr ! WorkflowManagerQueryFailure(e)
    }
  }

  private def replyToStatus(id: WorkflowId): Unit = {
    val sndr = sender()
    globalDataAccess.getWorkflowState(id) onComplete {
      // A "successful" return from the database API with a None value is actually a failure.  All legitimate
      // workflow IDs would have a status.
      case Success(None) => sndr ! WorkflowManagerStatusFailure(id, new WorkflowNotFoundException(s"Workflow $id not found"))
      case Success(s) => sndr ! WorkflowManagerStatusSuccess(id, s.get)
      case Failure(e) => sndr ! WorkflowManagerStatusFailure(id, e)
    }
  }

  startWith(Running, WorkflowManagerData(workflows = Map.empty))

  when (Running) {
    case Event(SubmitWorkflow(source), _) =>
      submitWorkflow(source, replyTo = Option(sender), id = None) map {
        case updatedEntry: WorkflowIdToActorRef => self ! AddEntryToWorkflowManagerData(updatedEntry)
      }
      //`recover` doesn't do anything since it's been handled already in the `submitWorkflow` function
      stay()
    case Event(WorkflowActorSubmitSuccess(replyTo, id), _) =>
      replyTo foreach { _ ! WorkflowManagerSubmitSuccess(id) }
      stay()
    case Event(WorkflowActorSubmitFailure(replyTo, e), _) =>
      replyTo foreach { _ ! WorkflowManagerSubmitFailure(e) }
      stay()
    case Event(WorkflowStatus(id), _) =>
      replyToStatus(id)
      stay()
    case Event(WorkflowQuery(rawParameters), _) =>
      replyToQuery(rawParameters)
      stay()
    case Event(WorkflowAbort(id), _) =>
      replyToAbort(id)
      stay()
    case Event(WorkflowOutputs(id), _) =>
      reply(id, workflowOutputs(id), WorkflowManagerWorkflowOutputsSuccess, WorkflowManagerWorkflowOutputsFailure)
      stay()
    case Event(CallOutputs(workflowId, callName), _) =>
      reply(workflowId, callName, callOutputs(workflowId, callName), WorkflowManagerCallOutputsSuccess, WorkflowManagerCallOutputsFailure)
      stay()
    case Event(CallStdoutStderr(workflowId, callName), _) =>
      val flatLogs = callStdoutStderr(workflowId, callName) map { _.flatten }
      reply(workflowId, callName, flatLogs, WorkflowManagerCallStdoutStderrSuccess, WorkflowManagerCallStdoutStderrFailure)
      stay()
    case Event(WorkflowStdoutStderr(id), _) =>
      val flatLogs = workflowStdoutStderr(id) map { _.mapValues(_.flatten) }
      reply(id, flatLogs, WorkflowManagerWorkflowStdoutStderrSuccess, WorkflowManagerWorkflowStdoutStderrFailure)
      stay()
    case Event(WorkflowMetadata(id), _) =>
      reply(id, workflowMetadata(id), WorkflowManagerWorkflowMetadataSuccess, WorkflowManagerWorkflowMetadataFailure)
      stay()
    case Event(SubscribeToWorkflow(id), data) =>
      //  NOTE: This fails silently. Currently we're ok w/ this, but you might not be in the future
      data.workflows.get(id) foreach {_ ! SubscribeTransitionCallBack(sender())}
      stay()
    case Event(RestartWorkflows(w :: ws), _) =>
      restartWorkflow(w) map {
        case updatedEntry: WorkflowIdToActorRef => self ! AddEntryToWorkflowManagerData(updatedEntry)
      }
      context.system.scheduler.scheduleOnce(RestartDelay) {
        self ! RestartWorkflows(ws)
      }
      stay()
    case Event(RestartWorkflows(Nil), _) =>
      // No more workflows need restarting.
      stay()
    case Event(CallCaching(id, parameters, callName), _) =>
      reply(id, callCaching(id, parameters, callName), WorkflowManagerCallCachingSuccess, WorkflowManagerCallCachingFailure)
      stay()
    case Event(Transition(workflowActor, _, toState: WorkflowState), data) if toState.isTerminal =>
      // Remove terminal actors from the store.
      val updatedData = data.remove(workflowActor)
      stay() using updatedData
    case Event(AbortAllWorkflows, data) if data.workflows.isEmpty =>
      goto(Done)
    case Event(AbortAllWorkflows, data) =>
      data.workflows.values.foreach { _ ! WorkflowActor.AbortWorkflow }
      goto(Aborting)
  }

  when (Aborting) {
    case Event(Transition(workflowActor, _, toState: WorkflowState), data) if toState.isTerminal =>
      // Remove this terminal actor from the workflowStore and log a progress message.
      val updatedData = data.remove(workflowActor)
      logger.info(s"$tag: Waiting for all workflows to abort (${updatedData.workflows.size} remaining).")
      // If there are no more workflows to abort we're done, otherwise just stay in the current state.
      (if (updatedData.workflows.isEmpty) goto(Done) else stay()) using updatedData
  }

  when (Done) { FSM.NullFunction }

  whenUnhandled {
    case Event(AddEntryToWorkflowManagerData(entry), _) =>
      logger.info(s"Updating WorkflowManager state. New Data: $entry")
      val newData = stateData.add(entry)
      stay() using newData
    // Uninteresting transition and current state notifications.
    case Event((Transition(_, _, _) | CurrentState(_, _)), _) => stay()
  }

  onTransition {
    case _ -> Done =>
      logger.info(s"$tag: All workflows aborted.")
      donePromise.trySuccess(())
  }

  /**
   * Returns a `Future[Any]` which will be failed if there is no workflow with the specified id.
   */
  private def assertWorkflowExistence(id: WorkflowId): Future[Any] = {
    // Confirm the workflow exists by querying its state.  If no state is found the workflow doesn't exist.
    globalDataAccess.getWorkflowState(id) map {
      case None => throw new WorkflowNotFoundException(s"Workflow '$id' not found")
      case _ =>
    }
  }

  private def assertCallExistence(id: WorkflowId, callFqn: FullyQualifiedName): Future[Any] = {
    globalDataAccess.getExecutionStatus(id, ExecutionDatabaseKey(callFqn, None, 1)) map {
      case None => throw new CallNotFoundException(s"Call '$callFqn' not found in workflow '$id'.")
      case _ =>
    }
  }

  private def workflowOutputs(id: WorkflowId): Future[engine.WorkflowOutputs] = {
    for {
      _ <- assertWorkflowExistence(id)
      outputs <- globalDataAccess.getWorkflowOutputs(id)
    } yield {
      SymbolStoreEntry.toWorkflowOutputs(outputs)
    }
  }

  private def callOutputs(workflowId: WorkflowId, callFqn: String): Future[core.CallOutputs] = {
    for {
      _ <- assertWorkflowExistence(workflowId)
      _ <- assertCallExistence(workflowId, callFqn)
      outputs <- globalDataAccess.getOutputs(workflowId, ExecutionDatabaseKey(callFqn, None, 1))
    } yield {
      SymbolStoreEntry.toCallOutputs(outputs)
    }
  }

  private def assertCallFqnWellFormed(descriptor: WorkflowDescriptor, callFqn: FullyQualifiedName): Future[Unit] = {
    descriptor.namespace.resolve(callFqn) match {
      case None => Future.failed(new UnsupportedOperationException(
        s"Expected a fully qualified name to have at least two parts but got $callFqn"))
      case _ => Future.successful(())
    }
  }

  private def callStdoutStderr(workflowId: WorkflowId, callFqn: String): Future[AttemptedCallLogs] = {
    for {
      _ <- assertWorkflowExistence(workflowId)
      descriptor <- globalDataAccess.getWorkflowExecutionAndAux(workflowId) flatMap workflowDescriptorFromExecutionAndAux
      _ <- assertCallExistence(workflowId, callFqn)
      _ <- assertCallFqnWellFormed(descriptor, callFqn)
      entries <- globalDataAccess.infosByExecution(workflowId, callFqn) map ExecutionInfosByExecution.toWorkflowLogs
    } yield entries.getOrElse(callFqn, List.empty)
  }

  private def workflowStdoutStderr(workflowId: WorkflowId): Future[WorkflowLogs] = {
    for {
      workflowState <- globalDataAccess.getWorkflowState(workflowId)
      // TODO: This assertion could also be added to the db layer: "In the future I'll fail if the workflow doesn't exist"
      _ <- WorkflowMetadataBuilder.assertWorkflowExistence(workflowId, workflowState)
      callLogOutputs <- globalDataAccess.infosByExecution(workflowId)
    } yield ExecutionInfosByExecution.toWorkflowLogs(callLogOutputs)
  }

  private def workflowMetadata(id: WorkflowId): Future[WorkflowMetadataResponse] = {
    WorkflowMetadataBuilder.workflowMetadata(id)
  }

  /** Submit the workflow and return an updated copy of the state data reflecting the addition of a
    * Workflow ID -> WorkflowActorRef entry.
    */
  private def submitWorkflow(source: WorkflowSourceFiles, replyTo: Option[ActorRef], id: Option[WorkflowId]): Future[WorkflowIdToActorRef] = {
    val workflowId: WorkflowId = id.getOrElse(WorkflowId.randomId())
    logger.info(s"$tag submitWorkflow input id = $id, effective id = $workflowId")
    val isRestart = id.isDefined

    import PromiseActor.EnhancedActorRef

    val message  = MaterializeWorkflowDescriptorActor.MaterializeWorkflow(workflowId, source, config)
    val materializeWorkflowDescriptorActor = context.actorOf(MaterializeWorkflowDescriptorActor.props(), name = s"MaterializeWorkflowDescriptorActor-$workflowId")
    materializeWorkflowDescriptorActor.askNoTimeout(message) map {
      case MaterializeWorkflowDescriptorSuccess(descriptor) =>
        val wfActor = context.actorOf(WorkflowActor.props(descriptor), name = s"WorkflowActor-$workflowId")
        wfActor ! SubscribeTransitionCallBack(self)
        wfActor ! (if (isRestart) Restart else Start(replyTo))
        logger.debug(s"Successfuly started ${wfActor.path} for Workflow $workflowId")
        context.stop(materializeWorkflowDescriptorActor)
        workflowId -> wfActor
      case MaterializeWorkflowDescriptorFailure(error) =>
        val messageOrBlank = Option(error.getMessage).mkString
        logger.error(error, s"$tag: Workflow failed submission: " + messageOrBlank)
        replyTo foreach { _ ! WorkflowManagerSubmitFailure(error)}
        // This error will be ignored currently as the only entity that needs to be notified about it
        // has already been sent a failure message above
        context.stop(materializeWorkflowDescriptorActor)
        throw error
    }
  }

  private def restartWorkflow(restartableWorkflow: WorkflowDescriptor): Future[WorkflowIdToActorRef] = {
    logger.info("Invoking restartableWorkflow on " + restartableWorkflow.id.shortString)
    submitWorkflow(restartableWorkflow.sourceFiles, replyTo = None, Option(restartableWorkflow.id))
  }

  private def restartIncompleteWorkflows(): Unit = {
    def logRestarts(restartableWorkflows: Traversable[WorkflowDescriptor]): Unit = {
      val num = restartableWorkflows.size
      val displayNum = if (num == 0) "no" else num.toString
      val plural = if (num == 1) "" else "s"

      logger.info(s"$tag Found $displayNum workflow$plural to restart.")

      if (num > 0) {
        val ids = restartableWorkflows.map { _.id.toString }.toSeq.sorted
        logger.info(s"$tag Restarting workflow ID$plural: " + ids.mkString(", "))
      }
    }

    val result = for {
      restartableWorkflowExecutionAndAuxes <- globalDataAccess.getWorkflowExecutionAndAuxByState(Seq(WorkflowSubmitted, WorkflowRunning))
      restartableWorkflows <- Future.sequence(restartableWorkflowExecutionAndAuxes map workflowDescriptorFromExecutionAndAux)
      _ = logRestarts(restartableWorkflows)
      _ = self ! RestartWorkflows(restartableWorkflows.toSeq)
    } yield ()

    result recover {
      case e: Throwable => logger.error(e, e.getMessage)
    }
  }

  private def query(rawParameters: Seq[(String, String)]): Future[WorkflowQueryResponse] = {
    for {
    // Future/Try to wrap the exception that might be thrown from WorkflowQueryParameters.apply.
      parameters <- Future.fromTry(Try(WorkflowQueryParameters(rawParameters)))
      response <- globalDataAccess.queryWorkflows(parameters)
    } yield response
  }

  private def callCaching(id: WorkflowId, parameters: QueryParameters, callName: Option[String]): Future[Int] = {
    for {
      _ <- assertWorkflowExistence(id)
      cachingParameters <- CallCachingParameters.from(id, callName, parameters)
      updateCount <- globalDataAccess.updateCallCaching(cachingParameters)
    } yield updateCount
  }
}
