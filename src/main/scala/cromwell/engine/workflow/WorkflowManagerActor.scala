package cromwell.engine.workflow

import akka.actor.FSM.{CurrentState, SubscribeTransitionCallBack, Transition}
import akka.actor._
import akka.event.Logging
import com.typesafe.config.ConfigFactory
import cromwell.engine
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.{Backend, CallLogs, CallMetadata}
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.db.slick._
import cromwell.engine.workflow.WorkflowActor.{Restart, Start}
import cromwell.engine.workflow.WorkflowManagerActor._
import cromwell.engine.{EnhancedFullyQualifiedName, _}
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice._
import lenthall.config.ScalaConfig.EnhancedScalaConfig
import org.joda.time.DateTime
import spray.json._
import wdl4s._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Future, Promise}
import scala.io.Source
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WorkflowManagerActor {
  class WorkflowNotFoundException(message: String) extends RuntimeException(message)
  class CallNotFoundException(message: String) extends RuntimeException(message)

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
  case object AbortAllWorkflows extends WorkflowManagerActorMessage

  def props(backend: Backend): Props = Props(new WorkflowManagerActor(backend))

  // FIXME hack to deal with one class of "singularity" where Cromwell isn't smart enough to launch only
  // as much work as can reasonably be handled.
  // How long to delay between restarting each workflow that needs to be restarted.  Attempting to
  // restart 500 workflows at exactly the same time crushes the database connection pool.
  lazy val RestartDelay = 200 milliseconds

  sealed trait WorkflowManagerState
  case object Running extends WorkflowManagerState
  case object Aborting extends WorkflowManagerState
  case object Done extends WorkflowManagerState

  type WorkflowActorRef = ActorRef

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
}

/**
 * Responses to messages:
 * SubmitWorkflow: Returns a `Future[WorkflowId]`
 * WorkflowStatus: Returns a `Future[Option[WorkflowState]]`
 * WorkflowOutputs: Returns a `Future[Option[binding.WorkflowOutputs]]` aka `Future[Option[Map[String, WdlValue]]`
 *
 */
class WorkflowManagerActor(backend: Backend) extends LoggingFSM[WorkflowManagerState, WorkflowManagerData] with CromwellActor {
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
      ConfigFactory.load.getConfig("backend").getBooleanOr("abortJobsOnTerminate", default = false)

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
      val updatedData = submitWorkflow(source, replyTo = Option(sender), id = None)
      stay() using updatedData
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
      val updatedData = restartWorkflow(w)
      context.system.scheduler.scheduleOnce(RestartDelay) {
        self ! RestartWorkflows(ws)
      }
      stay() using updatedData
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
    globalDataAccess.getExecutionStatus(id, ExecutionDatabaseKey(callFqn, None)) map {
      case None => throw new CallNotFoundException(s"Call '$callFqn' not found in workflow '$id'.")
      case _ =>
    }
  }

  /**
   * Retrieve the entries that produce stdout and stderr.
   */
  private def getCallLogKeys(id: WorkflowId, callFqn: FullyQualifiedName): Future[Seq[ExecutionDatabaseKey]] = {
    globalDataAccess.getExecutionStatuses(id, callFqn) map {
      case map if map.isEmpty => throw new CallNotFoundException(s"Call '$callFqn' not found in workflow '$id'.")
      case entries =>
        val callKeys = entries.keys filterNot { _.isCollector(entries.keys) }
        callKeys.toSeq.sortBy(_.index)
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

  private def callOutputs(workflowId: WorkflowId, callFqn: String): Future[engine.CallOutputs] = {
    for {
      _ <- assertWorkflowExistence(workflowId)
      _ <- assertCallExistence(workflowId, callFqn)
      outputs <- globalDataAccess.getOutputs(workflowId, ExecutionDatabaseKey(callFqn, None))
    } yield {
      SymbolStoreEntry.toCallOutputs(outputs)
    }
  }

  private def assertCallFqnWellFormed(descriptor: WorkflowDescriptor, callFqn: FullyQualifiedName): Try[String] = {
    descriptor.namespace.resolve(callFqn) match {
      case Some(c: Call) => Success(c.unqualifiedName)
      case _ => Failure(new UnsupportedOperationException(s"Expected a fully qualified name to have at least two parts but got $callFqn"))
    }
  }

  private def hasLogs(entries: Iterable[ExecutionDatabaseKey])(key: ExecutionDatabaseKey) = {
    !key.fqn.isScatter && !key.isCollector(entries)
  }

  private def callStdoutStderr(workflowId: WorkflowId, callFqn: String): Future[Seq[Seq[CallLogs]]] = {
    def callKey(descriptor: WorkflowDescriptor, callName: String, key: ExecutionDatabaseKey) =
      BackendCallKey(descriptor.namespace.workflow.findCallByName(callName).get, key.index, key.attempt)
    def backendCallFromKey(descriptor: WorkflowDescriptor, callName: String, key: ExecutionDatabaseKey) =
      backend.bindCall(descriptor, callKey(descriptor, callName, key))

    // Local import for FullyQualifiedName.isFinalCall since FullyQualifiedName is really String
    import cromwell.engine.finalcall.FinalCall._
    for {
        _ <- assertWorkflowExistence(workflowId)
        descriptor <- globalDataAccess.getWorkflow(workflowId)
        _ <- assertCallExistence(workflowId, callFqn)
        callName <- Future.fromTry(assertCallFqnWellFormed(descriptor, callFqn)) if !callFqn.isFinalCall
        callLogKeys <- getCallLogKeys(workflowId, callFqn)
        backendKeys <- Future.successful(callLogKeys.map(key => backendCallFromKey(descriptor, callName, key)))
      } yield backendKeys.groupBy(_.key.index).values.toSeq.sortBy(_.head.key.index) map { _.sortBy(_.key.attempt).map(_.stdoutStderr) }
  }

  private def workflowStdoutStderr(workflowId: WorkflowId): Future[Map[FullyQualifiedName, Seq[Seq[CallLogs]]]] = {
    def logMapFromStatusMap(descriptor: WorkflowDescriptor, statusMap: Map[ExecutionDatabaseKey, ExecutionStatus]): Try[Map[FullyQualifiedName, Seq[Seq[CallLogs]]]] = {
      // Local import for FullyQualifiedName.isFinalCall since FullyQualifiedName is really String
      import cromwell.engine.finalcall.FinalCall._
      Try {
        val sortedMap = statusMap.toSeq.sortBy(_._1.index)
        val callsToPaths = for {
          (key, status) <- sortedMap if !key.fqn.isFinalCall && hasLogs(statusMap.keys)(key)
          callName = assertCallFqnWellFormed(descriptor, key.fqn).get
          callKey = BackendCallKey(descriptor.namespace.workflow.findCallByName(callName).get, key.index, key.attempt)
          // TODO There should be an easier way than going as far as backend.bindCall just to retrieve stdout/err path
          backendCall = backend.bindCall(descriptor, callKey)
          callStandardOutput = backend.stdoutStderr(backendCall)
        } yield key -> callStandardOutput

        /* Some FP "magic" to transform the pairs of (key, logs) into the final result: grouped by FQNS, ordered by shards, and then ordered by attempts */
        callsToPaths groupBy { _._1.fqn } mapValues { key => key.groupBy(_._1.index).values.toSeq.sortBy(_.head._1.index) map { _.sortBy(_._1.attempt).map(_._2) }  }
      }
    }

    for {
      _ <- assertWorkflowExistence(workflowId)
      descriptor <- globalDataAccess.getWorkflow(workflowId)
      callToStatusMap <- globalDataAccess.getExecutionStatuses(workflowId)
      statusMap = callToStatusMap mapValues { _.executionStatus }
      callToLogsMap <- Future.fromTry(logMapFromStatusMap(descriptor, statusMap))
    } yield callToLogsMap
  }

  private def buildWorkflowMetadata(workflowExecution: WorkflowExecution,
                                    workflowExecutionAux: WorkflowExecutionAux,
                                    workflowOutputs: engine.WorkflowOutputs,
                                    callMetadata: Map[FullyQualifiedName, Seq[CallMetadata]]): WorkflowMetadataResponse = {

    val startDate = new DateTime(workflowExecution.startDt)
    val endDate = workflowExecution.endDt map { new DateTime(_) }
    val workflowInputs = Source.fromInputStream(workflowExecutionAux.jsonInputs.getAsciiStream).mkString.parseJson.asInstanceOf[JsObject]

    WorkflowMetadataResponse(
      id = workflowExecution.workflowExecutionUuid.toString,
      workflowName = workflowExecution.name,
      status = workflowExecution.status,
      // We currently do not make a distinction between the submission and start dates of a workflow, but it's
      // possible at least theoretically that a workflow might not begin to execute immediately upon submission.
      submission = startDate,
      start = Option(startDate),
      end = endDate,
      inputs = workflowInputs,
      outputs = Option(workflowOutputs) map { _.mapToValues },
      calls = callMetadata)
  }

  private def workflowMetadata(id: WorkflowId): Future[WorkflowMetadataResponse] = {
    for {
      workflowExecution <- globalDataAccess.getWorkflowExecution(id)
      workflowOutputs <- workflowOutputs(id)
      // The workflow has been persisted in the DB so we know the workflowExecutionId must be non-null,
      // so the .get on the Option is safe.
      workflowExecutionAux <- globalDataAccess.getWorkflowExecutionAux(WorkflowId.fromString(workflowExecution.workflowExecutionUuid))
      callStandardStreamsMap <- workflowStdoutStderr(id)
      callInputs <- globalDataAccess.getAllInputs(id)
      callOutputs <- globalDataAccess.getAllOutputs(id)
      infosByExecution <- globalDataAccess.infosByExecution(id)
      executionEvents <- globalDataAccess.getAllExecutionEvents(id)

      callMetadata = CallMetadataBuilder.build(infosByExecution, callStandardStreamsMap, callInputs, callOutputs, executionEvents)
      workflowMetadata = buildWorkflowMetadata(workflowExecution, workflowExecutionAux, workflowOutputs, callMetadata)

    } yield workflowMetadata
  }

  /** Submit the workflow and return an updated copy of the state data reflecting the addition of a
    * Workflow ID -> WorkflowActorRef entry.
    */
  private def submitWorkflow(source: WorkflowSourceFiles, replyTo: Option[ActorRef], id: Option[WorkflowId]): WorkflowManagerData = {
    val workflowId: WorkflowId = id.getOrElse(WorkflowId.randomId())
    logger.info(s"$tag submitWorkflow input id = $id, effective id = $workflowId")
    val isRestart = id.isDefined

    case class ActorAndStateData(actor: WorkflowActorRef, data: WorkflowManagerData)

    val actorAndData = for {
      descriptor <- Try(WorkflowDescriptor(workflowId, source))
      actor = context.actorOf(WorkflowActor.props(descriptor, backend), s"WorkflowActor-$workflowId")
      _ = actor ! SubscribeTransitionCallBack(self)
      data = stateData.add(workflowId -> actor)
    } yield ActorAndStateData(actor, data)

    actorAndData match {
      case Failure(e) =>
        val messageOrBlank = Option(e.getMessage).mkString
        logger.error(e, s"$tag: Workflow failed submission: " + messageOrBlank)
        replyTo foreach { _ ! WorkflowManagerSubmitFailure(e)}
        // Return the original state data if the workflow failed submission.
        stateData
      case Success(a) =>
        a.actor ! (if (isRestart) Restart else Start(replyTo))
        a.data
    }
  }

  private def restartWorkflow(restartableWorkflow: WorkflowDescriptor): WorkflowManagerData = {
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
      restartableWorkflows <- globalDataAccess.getWorkflowsByState(Seq(WorkflowSubmitted, WorkflowRunning))
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
