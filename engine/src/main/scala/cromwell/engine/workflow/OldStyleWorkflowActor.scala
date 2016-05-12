package cromwell.engine.workflow

import java.sql.SQLException

import akka.actor._
import akka.event.Logging
import akka.pattern.pipe
import cromwell.backend.{ExecutionEventEntry, ExecutionHash}
import cromwell.core.{CallOutput, CallOutputs}
import cromwell.database.obj.{Execution, ExecutionInfo}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus.{ExecutionStatus, _}
import cromwell.engine.backend._
import cromwell.engine.callactor.OldStyleCallActor
import cromwell.engine.callactor.OldStyleCallActor.CallActorMessage
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.EngineConverters.EnhancedExecution
import cromwell.engine.db.{CallStatus, ExecutionDatabaseKey, ExecutionInfosByExecution}
import cromwell.engine.finalcall.OldStyleFinalCall
import cromwell.engine.workflow.OldStyleWorkflowActor._
import cromwell.engine.workflow.OldStyleWorkflowManagerActor.{WorkflowActorSubmitFailure, WorkflowActorSubmitSuccess}
import cromwell.engine.{HostInputs, _}
import cromwell.logging.WorkflowLogger
import cromwell.util.TerminalUtil
import cromwell.webservice.WorkflowMetadataQueryParameters
import org.joda.time.DateTime
import wdl4s.types.WdlArrayType
import wdl4s.values.{WdlArray, WdlCallOutputsObject, WdlValue}
import wdl4s.{Scope, _}

import scala.annotation.tailrec
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleWorkflowActor {
  sealed trait WorkflowActorMessage
  case object GetFailureMessage extends WorkflowActorMessage
  case object AbortWorkflow extends WorkflowActorMessage
  sealed trait CallMessage extends WorkflowActorMessage {
    def callKey: ExecutionStoreKey
  }
  case class CallStarted(callKey: OutputKey, maybeCallLogs: Option[CallLogs]) extends CallMessage
  sealed trait TerminalCallMessage extends CallMessage
  case class CallAborted(callKey: OutputKey) extends TerminalCallMessage
  case class CallCompleted(callKey: OutputKey, callOutputs: CallOutputs, executionEvents: Seq[ExecutionEventEntry], returnCode: Int, hash: Option[ExecutionHash], resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor]) extends TerminalCallMessage
  case class ScatterCompleted(callKey: ScatterKey) extends TerminalCallMessage
  case class CallFailedRetryable(callKey: OutputKey, executionEvents: Seq[ExecutionEventEntry], returnCode: Option[Int], failure: Throwable) extends TerminalCallMessage
  case class CallFailedNonRetryable(callKey: OutputKey, executionEvents: Seq[ExecutionEventEntry], returnCode: Option[Int], failure: String) extends TerminalCallMessage
  case class CallFailedToInitialize(callKey: ExecutionStoreKey, reason: String) extends TerminalCallMessage
  case object Terminate extends WorkflowActorMessage
  final case class CachesCreated(startMode: StartMode) extends WorkflowActorMessage
  final case class AsyncFailure(t: Throwable) extends WorkflowActorMessage
  final case class PerformTransition(toState: WorkflowState) extends WorkflowActorMessage
  sealed trait PersistenceMessage extends WorkflowActorMessage {
    def callKey: ExecutionStoreKey
    def executionStatus: ExecutionStatus
  }

  final private case class PersistStatus(callKey: OutputKey, status: ExecutionStatus, callOutputs: Option[CallOutputs], returnCode: Option[Int],
                                         hash: Option[ExecutionHash], resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor], sender: ActorRef, message: TerminalCallMessage)
  final case class PersistenceSucceeded(callKey: ExecutionStoreKey,
                                        executionStatus: ExecutionStatus,
                                        outputs: Option[CallOutputs] = None) extends PersistenceMessage
  final case class PersistenceFailed(callKey: ExecutionStoreKey, executionStatus: ExecutionStatus) extends PersistenceMessage
  /** Used for exploded scatters which create many shards in one shot. */
  final case class PersistencesCompleted(callKeys: Traversable[ExecutionStoreKey],
                                         executionStatus: ExecutionStatus) extends WorkflowActorMessage
  /**
   * This message is sent from the context of the onTransition handler to trigger a `startRunnableCalls` invocation.
   * `startRunnableCalls` determines the calls which are pending persistence, information which is added to the
   * state data and passed forward to subsequent message processing.
   */
  case object StartRunnableCalls extends WorkflowActorMessage

  /**
    * Trigger a Running workflow to check whether or not it is complete (either successfully or not).
    */
  case object CheckForWorkflowComplete extends WorkflowActorMessage

  /**
   * Message sent to self to actually start a call.  This is necessary to synchronize access to the symbol cache
   * during call input retrieval.  The symbol cache is mutable state which should only be accessed in message
   * processing threads.
   * For initial start, assumes an execution is already persisted in Starting.  A restarted/resumed call should
   * be in Running, but the message handler will freshly persist the execution to Starting to side effect writing
   * a new call start time. */
  sealed trait CallStartMessage extends WorkflowActorMessage {
    def callKey: CallKey
    def startMode: OldStyleCallActor.StartMode
  }

  /** Represents starting a call for the first time, as opposed to a restart. */
  final case class InitialStartCall(override val callKey: CallKey,
                                    override val startMode: OldStyleCallActor.StartMode) extends CallStartMessage

  /** This signifies using an existing previously run call to fulfill the results of the callKey. */
  final case class UseCachedCall(override val callKey: BackendCallKey,
                                 override val startMode: OldStyleCallActor.UseCachedCall) extends CallStartMessage

  /** Represents restarting a call for backends which support restart. */
  final case class RestartCall(override val callKey: CallKey, override val startMode: OldStyleCallActor.StartMode) extends CallStartMessage

  sealed trait StartMode {
    def runInitialization(actor: OldStyleWorkflowActor): Future[Unit]
    def start(actor: OldStyleWorkflowActor): Unit
    def replyTo: Option[ActorRef]
  }

  implicit class EnhancedExecutionStoreKey(val key: ExecutionStoreKey) extends AnyVal {
    def toDatabaseKey: ExecutionDatabaseKey = ExecutionDatabaseKey(key.scope.fullyQualifiedName, key.index, key.attempt)
  }

  case class Start(replyTo: Option[ActorRef] = None) extends WorkflowActorMessage with StartMode {
    override def runInitialization(actor: OldStyleWorkflowActor): Future[Unit] = {
      // This only does the initialization for a newly created workflow.  For a restarted workflow we should be able
      // to assume the adjusted symbols already exist in the DB, but is it safe to assume the staged files are in place?
      actor.initializeWorkflow match {
        case Success(_) => actor.createWorkflow
        case Failure(ex) => Future.failed(ex)
      }
    }

    override def start(actor: OldStyleWorkflowActor) = actor.self ! StartRunnableCalls
  }

  case object Restart extends WorkflowActorMessage with StartMode {

    private def executionInfosToMap(infos: Seq[ExecutionInfo]): Map[String, Option[String]] = {
      infos.map(i => i.key -> i.value).toMap
    }

    private def restartOrResume(actor: OldStyleWorkflowActor, runningExecutions: Map[Execution, Seq[ExecutionInfo]]): Future[Unit] = {
      val executionsByCallFqn = runningExecutions.keys.groupBy(_.callFqn)

      def isRunningCollector(execution: Execution): Boolean = {
        executionsByCallFqn.get(execution.callFqn) match {
          case Some(xs) if xs.size > 1 && execution.index.toIndex.isEmpty && execution.executionStatus == ExecutionStatus.Running => true
          case _ => false
        }
      }

      val restartOrResumeAttempts = runningExecutions map { case (execution, executionInfos) =>
        val infos = executionInfosToMap(executionInfos)
        execution.toBackendCallKey(actor.workflow.namespace) match {
          case Success(backendCallKey) =>
            if (actor.backend.isResumable(backendCallKey, infos)) {
              Future.successful(())
            } else if (isRunningCollector(execution) || actor.backend.isRestartable(backendCallKey, infos)) {
              globalDataAccess.updateStatus(actor.workflow.id, Seq(execution.toKey), ExecutionStatus.NotStarted)
            } else {
              Future.failed(new Exception(s"Execution ${execution.executionId} is neither restartable or resumable"))
            }
          case Failure(ex) => Future.failed(ex)
        }
      }

      Future.sequence(restartOrResumeAttempts).map(_ => ())
    }

    private def filterResumableExecutions(actor: OldStyleWorkflowActor, runningExecutions: Map[Execution, Seq[ExecutionInfo]]): Future[Map[Execution, Seq[ExecutionInfo]]] = {
      val x = runningExecutions map { case (execution, executionInfos) =>
        val infos = executionInfosToMap(executionInfos)
        execution.toBackendCallKey(actor.workflow.namespace) match {
          case Success(backendCallKey) if actor.backend.isResumable(backendCallKey, infos) =>
            Future.successful(Option(execution -> executionInfos))
          case Success(_) => Future.successful(None)
          case Failure(ex) => Future.failed(ex)
        }
      }
      Future.sequence(x).map(_.flatten).map(_.toMap)
    }

    override def runInitialization(actor: OldStyleWorkflowActor): Future[Unit] = {
      def isWorkflowRestartable(allExecutions: Traversable[Execution]): Future[Unit] = {
        allExecutions.filter(e => Seq(ExecutionStatus.Aborted, ExecutionStatus.Failed).contains(e.executionStatus)) match {
          case es if es.nonEmpty =>
            val string = es.toSeq.sortWith((lt, rt) => lt.callFqn < rt.callFqn || (lt.callFqn == rt.callFqn && lt.index < rt.index)).mkString(" ")
            Future.failed(new Throwable(s"Workflow ${actor.workflow.id.shortString}: Cannot restart, found Failed and/or Aborted executions: $string"))
          case _ => Future.successful(())
        }
      }

      for {
        allExecutions <- globalDataAccess.getExecutions(actor.workflow.id)
        _ <- isWorkflowRestartable(allExecutions)
        runningExecutionsAndInfos <- globalDataAccess.runningExecutionsAndExecutionInfos(actor.workflow.id)
        executionsAndInfosMap = runningExecutionsAndInfos.map(x => x.execution -> x.executionInfos).toMap
        _ <- restartOrResume(actor, executionsAndInfosMap)
      } yield ()
    }

    override def start(actor: OldStyleWorkflowActor) = {
      val resumptionWork = for {
        runningExecutionsAndInfos <- globalDataAccess.runningExecutionsAndExecutionInfos(actor.workflow.id)
        executionsAndInfosAsMap = runningExecutionsAndInfos.map(x => x.execution -> x.executionInfos).toMap
        resumableExecutionsAndInfos <- filterResumableExecutions(actor, executionsAndInfosAsMap)
        _ = resumableExecutionsAndInfos foreach { case (exec, executionInfos) =>
          val executionInfosAsMap = executionInfosToMap(executionInfos)
          actor.self ! RestartCall(exec.toBackendCallKey(actor.workflow.namespace).get, OldStyleCallActor.Resume(executionInfosAsMap))
        }
        _ = actor.self ! StartRunnableCalls
      } yield ()

      resumptionWork onFailure {
        case t => actor.self ! AsyncFailure(t)
      }
    }

    override def replyTo: Option[ActorRef] = None
  }

  def props(descriptor: OldStyleWorkflowDescriptor): Props = {
    Props(OldStyleWorkflowActor(descriptor))
  }

  case class WorkflowData(startMode: Option[StartMode] = None,
                          pendingExecutions: Map[ExecutionStoreKey, Set[ExecutionStatus]] = Map.empty,
                          processingExecutions: Set[ExecutionStoreKey] = Set.empty) {

    def addPersisting(key: ExecutionStoreKey, status: ExecutionStatus)(implicit logger: WorkflowLogger) = {
      // Pending executions are modeled as a multi-valued map to deal with the case of newly created scatter shards.
      // For this particular case it's possible to have pending writes for both NotStarted and Starting states.
      // Adding message retries for executions in NotStarted state could significantly complicate startRunnableCalls
      // and make it less concurrent, and it's not particularly dire if those persists complete out of order.
      if (pendingExecutions.get(key).contains(status)) {
        logger.error(s"Unexpected add collision for key/status: ${key.tag}/$status")
      }
      val newPendingExecutions = pendingExecutions + (key -> (pendingExecutions.getOrElse(key, Set.empty) + status))
      this.copy(pendingExecutions = newPendingExecutions)
    }

    def removePersisting(key: ExecutionStoreKey, status: ExecutionStatus)(implicit logger: WorkflowLogger) = {
      if (!pendingExecutions.getOrElse(key, Set.empty).contains(status)) {
        logger.error(s"Unexpected remove of nonexistent key/status: ${key.tag}/$status")
      }
      val newStatusSet = pendingExecutions.getOrElse(key, Set.empty) - status
      // Remove the key/value completely if the status set is now empty.
      val newPendingExecutions = if (newStatusSet.isEmpty) pendingExecutions - key else pendingExecutions + (key -> newStatusSet)
      this.copy(pendingExecutions = newPendingExecutions)
    }

    def removePersisting(keys: Traversable[ExecutionStoreKey], status: ExecutionStatus)(implicit logger: WorkflowLogger): WorkflowData = {
      keys.foldLeft(this) { case (d, k) => d.removePersisting(k, status) }
    }

    def isPersistedRunning(key: ExecutionStoreKey)(implicit executionStore: ExecutionStore): Boolean = {
      executionStore.get(key).contains(ExecutionStatus.Running) && !pendingExecutions.contains(key)
    }

    def addProcessing(key: ExecutionStoreKey) = {
      val newProcessingExecutions = processingExecutions + key
      this.copy(processingExecutions = newProcessingExecutions)
    }

    def removeProcessing(key: ExecutionStoreKey)(implicit logger: WorkflowLogger): WorkflowData = {
      val newProcessingExecutions = processingExecutions - key
      this.copy(processingExecutions = newProcessingExecutions)
    }

    def isPending(key: ExecutionStoreKey) = pendingExecutions.contains(key)
    def isProcessing(key: ExecutionStoreKey) = processingExecutions.contains(key)
  }

  val AkkaTimeout = 5 seconds

  type ExecutionStore = Map[ExecutionStoreKey, ExecutionStatus]
  implicit class EnhancedStore(val es: ExecutionStore) extends AnyVal {
    // An entry is Done if one of its attempts is Done
    def isDone(e: ExecutionStoreEntry) = es exists { case (k, s) => k.scope == e._1.scope && k.index == e._1.index && s == ExecutionStatus.Done }
  }
  type ExecutionStoreEntry = (ExecutionStoreKey, ExecutionStatus)
  case class SymbolCacheKey(scopeName: String, input: Boolean)
  type SymbolCache = Map[SymbolCacheKey, Traversable[SymbolStoreEntry]]

  val TerminalStates = Vector(ExecutionStatus.Failed, ExecutionStatus.Done, ExecutionStatus.Aborted)

  def isExecutionStateFinished(es: ExecutionStatus): Boolean = TerminalStates contains es

  def isTerminal(status: ExecutionStatus): Boolean = TerminalStates contains status
  def isShard(key: BackendCallKey): Boolean = key.index.isDefined

  private val MarkdownMaxColumnChars = 100
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleWorkflowActor(workflow: OldStyleWorkflowDescriptor)
  extends LoggingFSM[WorkflowState, WorkflowData] with CromwellActor {

  implicit val actorSystem = context.system
  val backend = workflow.backend
  val workflowFailureMode = workflow.workflowFailureMode
  /*
   * This will not survive incoming Call-Scopification of backends but allows for use of engine functions everywhere for now,
   * in particular in the { input: ... } stanza of a call.
   */
  val engineFunctions = backend.engineFunctions(backend.fileSystems(workflow.workflowOptions), workflow.wfContext)

  def createWorkflow: Future[Unit] = {
    val symbolStoreEntries = buildSymbolStoreEntries(workflow, workflow.actualInputs)
    symbolCache = symbolStoreEntries.groupBy(entry => SymbolCacheKey(entry.scope, entry.isInput))
    val finalCalls = OldStyleFinalCall.createFinalCalls(workflow)
    globalDataAccess.createWorkflow(
      workflow, symbolStoreEntries, workflow.namespace.workflow.children ++ finalCalls, backend)
  }

  // This is passed as an implicit parameter to methods of classes in the companion object.
  // The execution and symbol caches are chunks of mutable state and must not be modified concurrently.
  implicit private var executionStore: ExecutionStore = _
  private var symbolCache: SymbolCache = _
  val akkaLogger = Logging(context.system, classOf[OldStyleWorkflowActor])
  implicit val logger: WorkflowLogger = WorkflowLogger("WorkflowActor", workflow, Option(akkaLogger))

  private[this] lazy val futureMetadata = new WorkflowMetadataBuilder(workflow.id,
    WorkflowMetadataQueryParameters(timings = false)).build()

  startWith(WorkflowSubmitted, WorkflowData())
  val startTime = System.nanoTime()
  /**
   * Try to generate output for a collector call, by collecting outputs for all of its shards.
   * It's fail-fast on shard output retrieval
   */
  private def generateCollectorOutput(collector: CollectorKey, shards: Iterable[BackendCallKey]): Try[CallOutputs] = Try {
    val shardsOutputs = shards.toSeq sortBy { _.index.fromIndex } map { e =>
      fetchCallOutputEntries(e) map { _.outputs } getOrElse(throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}"))
    }
    collector.scope.task.outputs map { taskOutput =>
      val wdlValues = shardsOutputs.map(s => s.getOrElse(taskOutput.name, throw new RuntimeException(s"Could not retrieve output ${taskOutput.name}")))
      val arrayOfValues = new WdlArray(WdlArrayType(taskOutput.wdlType), wdlValues)
      taskOutput.name -> CallOutput(arrayOfValues, workflow.hash(arrayOfValues))
    } toMap
  }

  private def findShardEntries(key: CollectorKey): Iterable[ExecutionStoreEntry] = executionStore collect {
    case (k: BackendCallKey, v) if k.scope == key.scope && isShard(k) => (k, v)
  }

  /** Messages self to perform the requested transition once all persistence is complete. */
  private def scheduleTransition(toState: WorkflowState): Unit = {
    // Called after the workflow Status is persisted, which means that if something fails here - Workflow will still be marked successful in the DB
    def handleTerminalWorkflow: Future[Unit] = {
      for {
          // FIXME If cleanUpForWorkflow fails - options are not updated and the terminate message is never sent
        _ <- backend.cleanUpForWorkflow(workflow)
        _ <- globalDataAccess.updateWorkflowOptions(workflow.id, workflow.workflowOptions.clearEncryptedValues)
        _ = self ! Terminate
      } yield ()
    }

    if (stateName != toState) {
      val transitionFuture = for {
      // Write the new workflow state before logging the change, tests assume the change is in effect when
      // the message is logged.
        _ <- globalDataAccess.updateWorkflowState(workflow.id, toState)
        _ = logger.info(s"Beginning transition from $stateName to $toState.")
        _ <- if (toState.isTerminal) handleTerminalWorkflow else Future.successful(())
      } yield ()

      transitionFuture recover {
        case e: Exception => logger.error(s"Failed to transition workflow status from $stateName to $toState", e)
      }

      // Only message self to perform the transition after the persistence of the new state is complete.
      transitionFuture map { _ => PerformTransition(toState) } pipeTo self
    }
  }

  /**
   * Captures the `ExecutionStoreKey`/`ExecutionStatus` combination for a call that has been created or started
   * as part of `startRunnableCalls`.  These are tracked to coordinate processing of messages based on completion
   * of persistence.
   */
  case class StartEntry(execution: ExecutionStoreKey, status: ExecutionStatus)

  /**
   * Collects the `StartEntry`s which have been persisted as part of starting runnable calls, as well as
   * any new entries which were added to the execution store (exploding scatters).
   */
  case class ExecutionStartResult(startEntries: Traversable[StartEntry],
                                  newEntries: Traversable[ExecutionStoreKey] = Seq.empty) {
    /** Merge two `ExecutionStartResult`s. */
    def plus(other: ExecutionStartResult): ExecutionStartResult = {
      ExecutionStartResult(startEntries = startEntries ++ other.startEntries, newEntries = newEntries ++ other.newEntries)
    }
  }

  private def initializeExecutionStore(startMode: StartMode): Future[Unit] = {
    val initializationCode = startMode.runInitialization(this)
    val futureCaches = for {
      _ <- initializationCode
      caches <- createCaches
    } yield caches

    futureCaches onComplete {
      case Success((executions, symbols)) =>
        executionStore = executions
        symbolCache = symbols
        self ! CachesCreated(startMode)
      case Failure(t) =>
        self ! AsyncFailure(t)
    }

    // "Lose" the actual value on purpose, the caller doesn't care/need to know about it
    futureCaches map { _ => () }
  }

  private def initializeWorkflow: Try[Unit] = backend.initializeForWorkflow(workflow)

  /**
   * Dump symbol and execution tables, start runnable calls, and message self to transition to the appropriate
   * next state.
   */
  private def dumpTables(): Future[Unit] = {
    for {
      symbols <- symbolsMarkdownTable
      _ = symbols foreach { table => logger.info(s"Initial symbols:\n\n$table") }
      executions <- executionsMarkdownTable
      _ = executions foreach { table => logger.info(s"Initial executions:\n\n$table") }
    } yield ()
  }

  when(WorkflowSubmitted) {
    case Event(startMode: StartMode, _) =>
      logger.info(s"$startMode message received")
      val sndr = sender()
      initializeExecutionStore(startMode) onComplete {
        case Success(id) => sndr ! WorkflowActorSubmitSuccess(startMode.replyTo, workflow.id)
        case Failure(e) => sndr ! WorkflowActorSubmitFailure(startMode.replyTo, e)
      }
      stay()
    case Event(CachesCreated(startMode), data) =>
      logger.info(s"ExecutionStoreCreated($startMode) message received")
      scheduleTransition(WorkflowRunning)
      stay() using data.copy(startMode = Option(startMode))
  }

  private def retryCall(callKey: ExecutionStoreKey,
                        data: WorkflowData) = {

    val keyClone = callKey.retryClone

    val insertCopy = globalDataAccess.insertCalls(workflow.id, List(keyClone), backend)

    insertCopy map { _ => PersistenceSucceeded(keyClone, ExecutionStatus.NotStarted) } recover {
      case e =>
        logger.error(s"Failed to clone ${callKey.tag} for retry attempt.", e)
        CheckForWorkflowComplete
    } pipeTo self

    data.addPersisting(keyClone, ExecutionStatus.NotStarted)
  }

  private def handleCallFailedRetryable(callKey: OutputKey,
                                        retryStatus: ExecutionStatus,
                                        returnCode: Option[Int],
                                        failure: Throwable,
                                        message: CallFailedRetryable,
                                        data: WorkflowData,
                                        events: Seq[ExecutionEventEntry]) = {
    val currentSender = sender()
    val persistEvents = globalDataAccess.setExecutionEvents(workflow.id, callKey.scope.fullyQualifiedName, callKey.index, callKey.attempt, events)

    persistEvents.failed.foreach {
      case e => logger.error(s"Failed to persist execution events for ${callKey.tag}.", e)
    }

    persistStatusThenAck(callKey, retryStatus, currentSender, message, None, returnCode, None, None)
    val updatedData = data.addPersisting(callKey, retryStatus)
    stay() using updatedData
  }


  private def handleCallCompleted(callKey: OutputKey,
                                  callOutputs: CallOutputs,
                                  executionEvents: Seq[ExecutionEventEntry],
                                  returnCode: Int,
                                  message: TerminalCallMessage,
                                  hash: Option[ExecutionHash],
                                  resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor],
                                  data: WorkflowData): State = {
    val currentSender = sender()
    val completionWork = for {
      // TODO These should be wrapped in a transaction so this happens atomically.
      _ <- globalDataAccess.setOutputs(workflow.id, callKey.toDatabaseKey, callOutputs,
        callKey.scope.rootWorkflow.outputs)
      _ <- globalDataAccess.setExecutionEvents(workflow.id, callKey.scope.fullyQualifiedName, callKey.index, callKey.attempt, executionEvents)
    } yield ()

    completionWork onComplete {
      case Failure(e) =>
        logger.error(s"Completion work failed for call ${callKey.tag}.", e)
        currentSender ! OldStyleCallActor.Ack(message)
        self ! CallFailedNonRetryable(callKey, executionEvents, Option(returnCode), e.getMessage)
      case Success(_) =>
        self ! PersistStatus(callKey, ExecutionStatus.Done, Option(callOutputs), Option(returnCode), hash, resultsClonedFrom, currentSender, message)
    }

    val updatedData = data.addProcessing(callKey)
    stay() using updatedData
  }

  /** Add the outputs for the specified `ExecutionStoreKey` to the symbol cache. */
  private def updateSymbolCache(executionKey: ExecutionStoreKey)(outputs: CallOutputs): Unit = {
    val newEntriesMap = outputs map { case (lqn, value) =>
      val storeKey = SymbolStoreKey(executionKey.scope.fullyQualifiedName, lqn, executionKey.index, input = false)
      new SymbolStoreEntry(storeKey, value.wdlValue.wdlType, Option(value.wdlValue), value.hash)
    } groupBy { entry => SymbolCacheKey(entry.scope, entry.isInput) }

    newEntriesMap foreach { case (key, entries) =>
      // SymbolCache is essentially a MultiMap, but that's a trait only for mutable Maps.
      symbolCache += key -> (symbolCache.getOrElse(key, Seq.empty) ++ entries)
    }
  }

  private def addCallFailureEvent(executionKey: ExecutionDatabaseKey, failureMessage: String): Unit = {
    logger.error(failureMessage)
    globalDataAccess.addCallFailureEvent(workflow.id, executionKey, FailureEventEntry(failureMessage, DateTime.now))
  }

  when(WorkflowRunning) {
    case Event(StartRunnableCalls, data) =>
      val updatedData = startRunnableCalls(data)
      self ! CheckForWorkflowComplete
      stay() using updatedData
    case Event(CallStarted(callKey, maybeCallLogs), data) if !data.isPending(callKey) =>
      executionStore += callKey -> ExecutionStatus.Running
      // TODO: We're not using/waiting for the future results
      // TODO: This should be in a single transaction!
      for {
        _ <- persistCallLogs(callKey, maybeCallLogs)
        _ <- persistStatus(callKey, ExecutionStatus.Running)
      } yield ()
      val updatedData = data.addPersisting(callKey, ExecutionStatus.Running)
      stay() using updatedData
    case Event(message: CallStarted, _) =>
      resendDueToPendingExecutionWrites(message)
      stay()
    case Event(message @ CallCompleted(callKey, outputs, executionEvents, returnCode, hash, resultsClonedFrom), data) if data.isPersistedRunning(callKey) && !data.isProcessing(callKey) =>
      handleCallCompleted(callKey, outputs, executionEvents, returnCode, message, hash, resultsClonedFrom, data)
    case Event(message @ CallCompleted(collectorKey: CollectorKey, outputs, executionEvents, returnCode, hash, resultsClonedFrom), data) if !data.isPending(collectorKey) =>
      // Collector keys are weird internal things and never go to Running state.
      handleCallCompleted(collectorKey, outputs, executionEvents, returnCode, message, hash, resultsClonedFrom, data)
    case Event(message @ CallCompleted(collectorKey: CollectorKey, _, _, _, _, _), _) =>
      resendDueToPendingExecutionWrites(message)
      stay()
    case Event(message @ PersistStatus(callKey, status, callOutputs, returnCode, hash, resultsClonedFrom, sender, terminalCallMessage), data) if !data.isPending(callKey) =>
      // This does not check the state of the FSM data, as this message can only be sent from the WA when and only when a Call status can be set to Success with 100% confidence
      // It also does not Ack to the sender as it can only be sent by the WA itself
      persistStatusThenAck(callKey, status, sender, terminalCallMessage, callOutputs, returnCode = returnCode, hash = hash, resultsClonedFrom = resultsClonedFrom)
      val updatedData = data.addPersisting(callKey, status)
      stay() using updatedData
    case Event(message @ CallFailedRetryable(callKey, events, returnCode, failure), data) if data.isPersistedRunning(callKey) && !data.isPending(callKey) =>
      addCallFailureEvent(callKey.toDatabaseKey, failure.getMessage)
      // Note: Currently Retryable == Preempted. If another reason than preemption arises that requires retrying at this level, this will need to be updated
      handleCallFailedRetryable(callKey, ExecutionStatus.Preempted, returnCode, failure, message, data, events)
    case Event(message @ CallFailedNonRetryable(callKey, events, returnCode, failure), data) if data.isPersistedRunning(callKey) =>
      persistStatusThenAck(callKey, ExecutionStatus.Failed, sender(), message, callOutputs = None, returnCode = returnCode)
      val updatedData = data.addPersisting(callKey, ExecutionStatus.Failed)
      addCallFailureEvent(callKey.toDatabaseKey, failure)
      stay() using updatedData
    case Event(message @ CallFailedToInitialize(callKey, reason), data) =>
      addCallFailureEvent(callKey.toDatabaseKey, s"Call failed to initialize: $reason")
      persistStatusThenAck(callKey, ExecutionStatus.Failed, sender(), message, callOutputs = None, returnCode = None)
      val updatedData = data.addPersisting(callKey, ExecutionStatus.Failed)
      stay() using updatedData
    case Event(message @ CallAborted(callKey), data) if data.isPersistedRunning(callKey) =>
      // Something funky's going on if aborts are coming through while the workflow's still running. But don't second-guess
      // by transitioning the whole workflow - the message is either still in the queue or this command was maybe
      // cancelled by some external system.
      executionStore += callKey -> ExecutionStatus.Aborted
      persistStatusThenAck(callKey, ExecutionStatus.Aborted, sender(), message)
      logger.warn(s"Call ${callKey.scope.unqualifiedName} was aborted but the workflow should still be running.")
      val updatedData = data.addPersisting(callKey, ExecutionStatus.Aborted)
      stay() using updatedData
    case Event(message: CallStartMessage, data) if !data.isPending(message.callKey) =>
      startCallWithMessage(message)
      stay()
    case Event(message: CallStartMessage, _) =>
      resendDueToPendingExecutionWrites(message)
      stay()
    case Event(PersistenceSucceeded(callKey, ExecutionStatus.Done, callOutputs), data) =>
      executionStore += callKey -> ExecutionStatus.Done
      callOutputs foreach updateSymbolCache(callKey)
      logger.debug(s"In state WorkflowRunning: Got PersistenceCompleted message for Done call ${callKey.tag}")
      val updatedData = startRunnableCalls(data)
      val finalData = updatedData.removePersisting(callKey, ExecutionStatus.Done).removeProcessing(callKey)
      self ! CheckForWorkflowComplete
      stay using finalData
    case Event(PersistenceSucceeded(callKey, ExecutionStatus.Preempted, callOutputs), data) =>
      executionStore += callKey -> ExecutionStatus.Preempted
      val updatedData = retryCall(callKey, data).removePersisting(callKey, ExecutionStatus.Preempted).removeProcessing(callKey)
      stay using updatedData
    case Event(PersistenceSucceeded(callKey, ExecutionStatus.Failed, callOutputs), data) =>
      executionStore += callKey -> ExecutionStatus.Failed
      val finalData = data.removePersisting(callKey, ExecutionStatus.Failed).removeProcessing(callKey)
      self ! CheckForWorkflowComplete
      stay using finalData
    case Event(PersistenceSucceeded(callKey, ExecutionStatus.NotStarted, callOutputs), data) =>
      executionStore += callKey -> ExecutionStatus.NotStarted
      val updatedData = startRunnableCalls(data).removePersisting(callKey, ExecutionStatus.NotStarted).removeProcessing(callKey)
      self ! CheckForWorkflowComplete
      stay using updatedData
    case Event(CheckForWorkflowComplete, data) =>
      checkForWorkflowComplete(data)
      stay()
  }

  private def checkForWorkflowComplete(data: WorkflowData): Unit = {
    val executionStatuses = executionStore.values.toSeq

    if (!executionStatuses.contains(ExecutionStatus.Running) && !executionStatuses.contains(ExecutionStatus.Starting) && data.pendingExecutions.isEmpty) {
      if (executionStatuses.contains(ExecutionStatus.Failed)) {
        scheduleTransition(WorkflowFailed)
      } else if (isWorkflowDone) {
        scheduleTransition(WorkflowSucceeded)
      }
    }
  }

  private def startBackendCallWithMessage(message: CallStartMessage, backendCallKey: BackendCallKey, callInputs: Map[String, WdlValue]) = {
      // TODO: The block of code is only run on non-shards because of issues with DSDEEPB-2490
      if (!backendCallKey.index.isShard) {
        val updateDbCallInputs = globalDataAccess.updateCallInputs(workflow.id, backendCallKey, callInputs)
        updateDbCallInputs onComplete {
          case Success(i) =>
            logger.debug(s"$i call input expression(s) updated in database.")
            startActor(backendCallKey, callInputs, message.startMode)
          case Failure(e) => self ! AsyncFailure(new SQLException(s"Failed to update symbol inputs for ${backendCallKey.scope.fullyQualifiedName}.${backendCallKey.tag}.${backendCallKey.index}", e))
        }
      } else {
        startActor(backendCallKey, callInputs, message.startMode)
      }
  }

  private def startCallWithMessage(message: CallStartMessage) = {
    val callKey = message.callKey
    callKey match {
      case backendCallKey: BackendCallKey => fetchLocallyQualifiedInputs(backendCallKey) match {
        case Success(callInputs) => startBackendCallWithMessage(message, backendCallKey, callInputs)
        case Failure(t) =>
          self ! CallFailedToInitialize(callKey, s"Failed to fetch locally qualified inputs: ${t.getMessage}")
      }
      case finalCallKey: FinalCallKey => startActor(finalCallKey, Map.empty, message.startMode)
    }
  }

  when(WorkflowSucceeded) { FSM.NullFunction }

  when(WorkflowFailed) { FSM.NullFunction }

  when(WorkflowAborted) { FSM.NullFunction }

  onTransition {
    case WorkflowSubmitted -> WorkflowRunning =>
      stateData.startMode.get.start(this)
  }

  when(WorkflowAborting) {
    case Event(message @ CallCompleted(callKey, outputs, executionEvents, returnCode, hash, resultsClonedFrom), data) if data.isPersistedRunning(callKey) =>
      handleCallCompleted(callKey, outputs, executionEvents, returnCode, message, hash, resultsClonedFrom, data)
    case Event(message @ CallCompleted(collectorKey: CollectorKey, outputs, executionEvents, returnCode, hash, resultsClonedFrom), data) if !data.isPending(collectorKey) =>
      // Collector keys are weird internal things and never go to Running state.
      handleCallCompleted(collectorKey, outputs, executionEvents, returnCode, message, hash, resultsClonedFrom, data)
    case Event(message @ CallAborted(callKey), data) if data.isPersistedRunning(callKey) =>
      executionStore += callKey -> ExecutionStatus.Aborted
      persistStatusThenAck(callKey, ExecutionStatus.Aborted, sender(), message)
      val updatedData = data.addPersisting(callKey, ExecutionStatus.Aborted)
      if (isWorkflowAborted) scheduleTransition(WorkflowAborted)
      stay() using updatedData
    case Event(message @ CallFailedRetryable(callKey, events, returnCode, failure), data) if data.isPersistedRunning(callKey) =>
      persistStatusThenAck(callKey, ExecutionStatus.Failed, sender(), message, callOutputs = None, returnCode = returnCode)
      val updatedData = data.addPersisting(callKey, ExecutionStatus.Failed)
      if (isWorkflowAborted) scheduleTransition(WorkflowAborted)
      stay() using updatedData
    case Event(PersistenceSucceeded(callKey, ExecutionStatus.Done, callOutputs), data) =>
      executionStore += callKey -> ExecutionStatus.Done
      callOutputs foreach updateSymbolCache(callKey)
      logger.debug(s"In state WorkflowAborting: Got PersistenceCompleted message for Done call ${callKey.tag}")
      if (isWorkflowAborted) scheduleTransition(WorkflowAborted)
      val updatedData = data.removePersisting(callKey, ExecutionStatus.Done)
      stay() using updatedData
  }

  /**
   * It is legitimate to switch states if the current state is not terminal and the target state is not the same as
   * the current state.
   */
  private def canSwitchTo(toState: WorkflowState): Boolean = {
    !stateName.isTerminal && stateName != toState
  }

  private def resendDueToPendingExecutionWrites(message: Any): Unit = {
    val ResendInterval = 100.milliseconds
    val numPersists = stateData.pendingExecutions.size
    val stringPersists = stateData.pendingExecutions map { case (key, statuses) => key.tag + "/" + statuses.mkString(", ") } mkString ";   "
    logger.debug(s"Rescheduling message $message with $ResendInterval delay due to $numPersists pending persists: $stringPersists")
    context.system.scheduler.scheduleOnce(ResendInterval, self, message)
  }

  private def removePendingCallKeyPersistence(data: WorkflowData, message: PersistenceMessage): FSM.State[WorkflowState, WorkflowData] = {
    val updatedData = data.removePersisting(message.callKey, message.executionStatus)
    self ! CheckForWorkflowComplete
    stay using updatedData
  }

  whenUnhandled {
    case Event(AbortWorkflow, _) =>
      context.children foreach { _ ! OldStyleCallActor.AbortCall }
      scheduleTransition(WorkflowAborting)
      stay()
    case Event(ScatterCompleted(scatterKey), data) =>
      // This case is common to WorkflowRunning and WorkflowAborting.
      log.debug(s"Got ScatterCompleted($scatterKey)")
      executionStore += scatterKey -> ExecutionStatus.Done
      persistStatus(scatterKey, ExecutionStatus.Done, callOutputs = None, returnCode = Option(0))
      val updatedData = data.addPersisting(scatterKey, ExecutionStatus.Done)
      stay() using updatedData
    case Event(callMessage: CallMessage, _) =>
      logger.debug(s"Dropping message for ineligible key: ${callMessage.callKey.tag}")
      // Running and Aborting have explicit handlers for Running and no pending writes, so either we are not in those
      // states or the specified key is no longer Running and/or has pending writes.
      stay()
    case Event(m: PersistenceMessage, data) =>
      logger.debug(s"Got whenUnhandled message: $m")
      removePendingCallKeyPersistence(data, m)
    case Event(PersistencesCompleted(callKeys, status), data) =>
      val stringKeys = callKeys map { _.tag } mkString ", "
      logger.debug(s"Got whenUnhandled PersistencesCompleted($stringKeys)")
      val updatedData = data.removePersisting(callKeys, status)
      stay using updatedData
    case Event(PerformTransition(toState), data) if data.pendingExecutions.isEmpty && canSwitchTo(toState) =>
      logger.info(s"transitioning from $stateName to $toState.")
      goto(toState)
    case Event(PerformTransition(toState), _) if !canSwitchTo(toState) =>
      // Drop the message.  This is a race condition with a retried PerformTransition message when
      // the FSM has transitioned to a terminal state between the original and retried messages.
      logger.debug(s"Dropping PerformTransition message: $stateName -> $toState")
      stay()
    case Event(message: PerformTransition, _) =>
      resendDueToPendingExecutionWrites(message)
      stay()
    case Event(AsyncFailure(t), data) if data.pendingExecutions.isEmpty =>
      logger.error(t.getMessage, t)
      self ! CheckForWorkflowComplete
      stay()
    case Event(message @ AsyncFailure(t), _) =>
      // This is the unusual combination of warn + throwable logging since the expectation is that this will eventually
      // be logged in the case above as an error, but if for some weird reason this actor never ends up in that
      // state we don't want to be completely blind to the cause of the AsyncFailure.
      logger.warn(t.getMessage, t)
      resendDueToPendingExecutionWrites(message)
      stay()
    case Event(Terminate, data) if data.pendingExecutions.isEmpty && stateName.isTerminal =>
      shutDown()
      stay()
    case Event(Terminate, _) =>
      // Don't actually terminate as there are pending writes and/or this is not yet transitioned to a terminal state.
      // Resend the message after a delay.
      resendDueToPendingExecutionWrites(Terminate)
      stay()
    case Event(e, _) =>
      logger.debug(s"received unhandled event $e while in state $stateName")
      stay()
  }

  private def persistStatus(storeKey: ExecutionStoreKey,
                            executionStatus: ExecutionStatus,
                            callOutputs: Option[CallOutputs] = None,
                            returnCode: Option[Int] = None,
                            hash: Option[ExecutionHash] = None,
                            resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor] = None): Future[Unit] = {

    logger.info(s"persisting status of ${storeKey.tag} to $executionStatus.")

    val persistFuture = executionStatus match {
      case ExecutionStatus.Starting => globalDataAccess.setStartingStatus(workflow.id, List(storeKey.toDatabaseKey))
      case terminal if terminal.isTerminal => globalDataAccess.setTerminalStatus(workflow.id, storeKey.toDatabaseKey, terminal, returnCode, hash, resultsClonedFrom)
      case other => globalDataAccess.updateStatus(workflow.id, List(storeKey.toDatabaseKey), other)
    }

    persistFuture onComplete {
      case Success(_) => self ! PersistenceSucceeded(storeKey, executionStatus, callOutputs)
      case Failure(t) =>
        logger.error(s"Error persisting status of call ${storeKey.tag} to $executionStatus", t)
        self ! CheckForWorkflowComplete
    }
    persistFuture
  }

  private def persistStatusThenAck(storeKey: ExecutionStoreKey,
                                   executionStatus: ExecutionStatus,
                                   recipient: ActorRef,
                                   message: TerminalCallMessage,
                                   callOutputs: Option[CallOutputs] = None,
                                   returnCode: Option[Int] = None,
                                   hash: Option[ExecutionHash] = None,
                                   resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor] = None): Unit = {

    persistStatus(storeKey, executionStatus, callOutputs, returnCode, hash, resultsClonedFrom) map { _ => OldStyleCallActor.Ack(message) } pipeTo recipient
  }

  private def persistCallLogs(storeKey: ExecutionStoreKey, maybeCallLogs: Option[CallLogs]): Future[Unit] = {
    (maybeCallLogs, storeKey) match {
      case (Some(callLogs), backendCallKey: BackendCallKey) =>
        globalDataAccess.upsertExecutionInfo(
          workflow.id,
          backendCallKey,
          ExecutionInfosByExecution.toCallLogMap(callLogs),
          actorSystem
        )
      case _ => Future.successful(Unit)
    }
  }

  private def startActor(callKey: CallKey, locallyQualifiedInputs: CallInputs,
                         callActorMessage: CallActorMessage): Unit = {
    if (locallyQualifiedInputs.nonEmpty) {
      val inputs = locallyQualifiedInputs map { case(lqn, value) => s"  $lqn -> $value" } mkString "\n"
      logger.info(s"inputs for call '${callKey.tag}':\n$inputs")
    } else {
      logger.info(s"no inputs for call '${callKey.tag}'")
    }

    val callActorName = s"CallActor-${workflow.id}-${callKey.tag}"
    val futureCallActorProps = callKey match {
      case backendCallKey: BackendCallKey =>
        // PBE This awful `workflow.copy` bit was inspired by the RetryableCallsSpec that passes a tweaked version
        // of the local backend down through WMA.  But it turns out that WorkflowDescriptor.apply doesn't use WMA's
        // value of backend, it talks directly to CromwellBackend and gets the regular local backend, which breaks the test.

        // This backend parameter to WMA is confusing and in the world of pluggable backends somehow makes even less
        // sense than using backends on individual workflows.  When we get real pluggable backends that test can
        // just set its special backend cleanly in config, and all this awfulness can go away.
        val descriptor = OldStyleBackendCallJobDescriptor(workflow.copy(backend = backend), backendCallKey, locallyQualifiedInputs)
        Future.successful(OldStyleCallActor.props(descriptor))
      case finalCallKey: FinalCallKey =>
        futureMetadata map { metadata => OldStyleCallActor.props(FinalCallJobDescriptor(workflow, finalCallKey, metadata)) }
    }

    futureCallActorProps onComplete {
      case Success(callActorProps) =>
        val callActor = context.actorOf(callActorProps, callActorName)
        callActor ! callActorMessage
        logger.info(s"created call actor for ${callKey.tag}.")
      case Failure(throwable) =>
        logger.error(s"failed to create call actor for ${callKey.tag}.", throwable)
        self ! CallFailedToInitialize(callKey,
          s"failed to create call actor for ${callKey.tag}: ${throwable.getMessage}")
    }
  }

  /**
    * Attempt to start all runnable calls and return updated state data.  This will create a new copy
    * of the state data including new pending persists.
    */
  @tailrec
  private def startRunnableCalls(data: WorkflowData): WorkflowData = {

    def upstreamEntries(entry: ExecutionStoreKey, prerequisiteScope: Scope): Seq[ExecutionStoreEntry] = {
      prerequisiteScope.closestCommonAncestor(entry.scope) match {
        /**
         * If this entry refers to a Scope which has a common ancestor with prerequisiteScope
         * and that common ancestor is a Scatter block, then find the shard with the same index
         * as 'entry'.  In other words, if you're in the same scatter block as your pre-requisite
         * scope, then depend on the shard (with same index).
         *
         * NOTE: this algorithm was designed for ONE-LEVEL of scattering and probably does not
         * work as-is for nested scatter blocks
         */
        case Some(ancestor: Scatter) =>
          executionStore filter { case(k, _) => k.scope == prerequisiteScope && k.index == entry.index } toSeq

        /**
         * Otherwise, simply refer to the entry the collector entry.  This means that 'entry' depends
         * on every shard of the pre-requisite scope to finish.
         */
        case _ =>
          executionStore filter { case(k, _) => k.scope == prerequisiteScope && k.index.isEmpty } toSeq
      }
    }

    def arePrerequisitesDone(key: ExecutionStoreKey): Boolean = {
      val upstream = key.scope.prerequisiteScopes.map(s => upstreamEntries(key, s))
      val downstream = key match {
        case collector: CollectorKey => findShardEntries(collector)
        case _ => Nil
      }
      val dependencies = upstream.flatten ++ downstream
      val dependenciesResolved = executionStore.filter(dependencies) forall executionStore.isDone

      /**
       * We need to make sure that all prerequisiteScopes have been resolved to some entry before going forward.
       * If a scope cannot be resolved it may be because it is in a scatter that has not been populated yet,
       * therefore there is no entry in the executionStore for this scope.
       * If that's the case this prerequisiteScope has not been run yet, hence the (upstream forall {_.nonEmpty})
       */
      (upstream forall { _.nonEmpty }) && dependenciesResolved
    }

    lazy val nonFinalCallsComplete = executionStore forall {
      case (key: FinalCallKey, _) => true
      case (_, x: ExecutionStatus) if x.isTerminal => true
      case _ => false
    }

    lazy val hasAnythingFailed = executionStore.values exists { x => x == Failed || x == Aborted }
    lazy val allowNewCalls = workflowFailureMode.allowNewCallsAfterFailure || !hasAnythingFailed

    def isRunnable(entry: ExecutionStoreEntry) = {
      entry match {
        case (key: FinalCallKey, ExecutionStatus.NotStarted) => nonFinalCallsComplete && allowNewCalls
        case (key, ExecutionStatus.NotStarted) => arePrerequisitesDone(key) && allowNewCalls
        case _ => false
      }
    }

    val runnableEntries = executionStore filter isRunnable

    val runnableCalls = runnableEntries collect { case(k: BackendCallKey, v) => k.scope }
    if (runnableCalls.nonEmpty)
      logger.info(s"starting calls: " + runnableCalls.map(_.fullyQualifiedName).toSeq.sorted.mkString(", "))

    val entries: Map[ExecutionStoreKey, Try[ExecutionStartResult]] = runnableEntries map {
      case (k: BackendCallKey, _) => k -> processRunnableCall(k)
      case (k: FinalCallKey, _) => k -> processRunnableCall(k)
      case (k: ScatterKey, _) => k -> processRunnableScatter(k)
      case (k: CollectorKey, _) => k -> processRunnableCollector(k)
      case (k, v) =>
        val message = s"Unknown entry in execution store:\nKEY: ${k.tag}\nVALUE:$v"
        logger.error(message)
        k -> Failure(new UnsupportedOperationException(message))
    }

    val zero = ExecutionStartResult(startEntries = Set.empty)

    // Fold across the entries. Successes get added to the accumulator. Failures cause a message.
    val result = entries.foldLeft(zero) {
      case (acc, (_, Success(x))) =>
        acc plus x
      case (acc, (key, Failure(e))) =>
        self ! CallFailedToInitialize(key, s"Failed to start call: ${e.getMessage}")
        acc
    }

    val updatedData = result.startEntries.foldLeft(data) { case (d, StartEntry(k, s)) => d.addPersisting(k, s) }
    if (result.newEntries.nonEmpty) startRunnableCalls(updatedData) else updatedData
  }

  private def lookupNamespace(name: String): Try[WdlNamespace] = {
    workflow.namespace.namespaces find { _.importedAs.contains(name) } match {
      case Some(x) => Success(x)
      case _ => Failure(new WdlExpressionException(s"Could not resolve $name as a namespace"))
    }
  }

  private def lookupCall(key: ExecutionStoreKey, workflow: Workflow)(name: String): Try[WdlCallOutputsObject] = {
    workflow.calls find { _.unqualifiedName == name } match {
      case Some(matchedCall) =>
        /**
         * After matching the Call, this determines if the `key` depends on a single shard
         * of a scatter'd job or if it depends on the whole thing.  Right now, the heuristic
         * is "If we're both in a scatter block together, then I depend on a shard.  If not,
         * I depend on the collected value"
         *
         * TODO: nested-scatter - this will likely not be sufficient for nested scatters
         */
        val index: ExecutionIndex = matchedCall.closestCommonAncestor(key.scope) flatMap {
          case s: Scatter => key.index
          case _ => None
        }
        fetchCallOutputEntries(findCallKey(matchedCall, index) getOrElse {
          throw new WdlExpressionException(s"Could not find a callKey for name '${matchedCall.unqualifiedName}'")
        })
      case None => Failure(new WdlExpressionException(s"Could not find a call with name '$name'"))
    }
  }

  private def lookupDeclaration(workflow: Workflow)(name: String): Try[WdlValue] = {
    workflow.scopedDeclarations find { _.name == name } match {
      case Some(declaration) => fetchFullyQualifiedName(declaration.fullyQualifiedName)
      case None => Failure(new WdlExpressionException(s"Could not find a declaration with name '$name'"))
    }
  }

  private def lookupScatterVariable(callKey: BackendCallKey, workflow: Workflow)(name: String): Try[WdlValue] = {
    val scatterBlock = callKey.scope.ancestry collect { case s: Scatter => s } find { _.item == name }
    val scatterCollection = scatterBlock map { s =>
      s.collection.evaluate(scatterCollectionLookupFunction(workflow, callKey), engineFunctions) match {
        case Success(v: WdlArray) if callKey.index.isDefined =>
          if (v.value.isDefinedAt(callKey.index.get))
            Success(v.value(callKey.index.get))
          else
            Failure(new WdlExpressionException(s"Index ${callKey.index.get} out of bounds for $name array."))
        case Success(v: WdlArray) => Failure(new WdlExpressionException(s"$name evaluated to an Array but $callKey has no index"))
        case _ => Failure(new WdlExpressionException(s"$name did not evaluate to a WdlArray"))
      }
    }
    scatterCollection.getOrElse(
      Failure(new WdlExpressionException(s"$name is does not reference a scattered variable"))
    )
  }

  private def resolveIdentifierOrElse(identifierString: String, resolvers: ((String) => Try[WdlValue]) *)(orElse: => Try[WdlValue]): WdlValue = {
    /* Try each of the resolver functions in order.  This uses a lazy Stream to only call a resolver function if a
     * preceding resolver function failed to resolve the identifier. */
    val attemptedResolutions = Stream(resolvers: _*) map { _(identifierString) } find { _.isSuccess }

    /* Return the first successful function's value or throw an exception. */
    attemptedResolutions.getOrElse(orElse).get
  }

  private def findCallKey(call: Call, index: ExecutionIndex): Option[OutputKey] = {
    executionStore.collect({
      case (k: OutputKey, _) if k.scope == call && k.index == index => k
    }).headOption
  }

  def fetchLocallyQualifiedInputs(callKey: BackendCallKey): Try[Map[String, WdlValue]] = Try {
    val parentWorkflow = callKey.scope.rootWorkflow

    def lookup(identifier: String): WdlValue = {
      /* This algorithm defines three ways to lookup an identifier in order of their precedence:
       *
       *   1) Traverse up the scope hierarchy and see if the variable reference any scatter item
       *   2) Look for a WdlNamespace with matching name
       *   3) Look for a Call with a matching name (perhaps using a scope resolution algorithm)
       *   4) Look for a Declaration with a matching name (perhaps using a scope resolution algorithm)
       *
       *  Each method is tried individually and the first to return a Success value takes precedence.
       */

      resolveIdentifierOrElse(identifier, lookupScatterVariable(callKey, parentWorkflow), lookupNamespace, lookupCall(callKey, parentWorkflow), lookupDeclaration(parentWorkflow)) {
        Failure(new WdlExpressionException(s"Could not resolve $identifier as a scatter variable, namespace, call, or declaration"))
      }
    }

    val entries = fetchCallInputEntries(callKey)
    entries.map { entry =>
      // .get are used below because the exception will be captured by the Try
      val declaration = findTaskDeclaration(entry.scope, entry.key.name).get

      val value = entry.wdlValue match {
        case Some(e: WdlExpression) => e.evaluate(lookup, engineFunctions)
        case Some(v) => Success(v)
        case _ => Failure(new WdlExpressionException("Unknown error"))
      }
      val coercedValue = value.flatMap(x => declaration.wdlType.coerceRawValue(x))
      entry.key.name -> coercedValue.get
    }.toMap
  }

  private def findTaskDeclaration(callFqn: String, inputName: String): Try[Declaration] = {
    val exception = new WdlExpressionException(s"Could not find task input '$inputName' for call '$callFqn'")
    workflow.namespace.resolve(callFqn) match {
      case Some(c:Call) =>
        c.task.declarations.find(_.name == inputName) match {
          case Some(decl) => Success(decl)
          case None => Failure(exception)
        }
      case _ => Failure(exception)
    }
  }

  /** Return all input symbols with this FQN. */
  private def inputSymbolsByFullyQualifiedName(fqn: FullyQualifiedName): Traversable[SymbolStoreEntry] = {
    val (scope, variable) = fqn.scopeAndVariableName
    symbolCache.getOrElse(SymbolCacheKey(scope, input = true), Seq.empty) filter { _.key.name == variable }
  }

  private def fetchFullyQualifiedName(fqn: FullyQualifiedName): Try[WdlValue] = {
    inputSymbolsByFullyQualifiedName(fqn) match {
      case t: Traversable[SymbolStoreEntry] if t.isEmpty =>
        Failure(new WdlExpressionException(s"Could not find a declaration with fully-qualified name '$fqn'"))
      case t: Traversable[SymbolStoreEntry] if t.size > 1 =>
        Failure(new WdlExpressionException(s"Expected only one declaration for fully-qualified name '$fqn', got ${t.size}"))
      case t: Traversable[SymbolStoreEntry] => t.head.wdlValue match {
        case Some(value) => Success(value)
        case None => Failure(new WdlExpressionException(s"No value defined for fully-qualified name $fqn"))
      }
    }
  }

  private def callOutputEntries(outputKey: OutputKey): Traversable[SymbolStoreEntry] = {
    val cacheKey = SymbolCacheKey(outputKey.scope.fullyQualifiedName, input = false)
    symbolCache.getOrElse(cacheKey, Seq.empty) filter { _.key.index == outputKey.index }
  }

  private def fetchCallOutputEntries(outputKey: OutputKey): Try[WdlCallOutputsObject] = {
    val callOutputsAsMap = callOutputEntries(outputKey).map(entry => entry.key.name -> entry.wdlValue).toMap
    callOutputsAsMap find { case (k, v) => v.isEmpty } match {
      case Some(noneValue) => Failure(new WdlExpressionException(s"Could not evaluate call ${outputKey.scope.unqualifiedName} because some of its inputs are not defined (i.e. ${noneValue._1}"))
      // TODO: .asInstanceOf[Call]?
      case None => Success(WdlCallOutputsObject(outputKey.scope.asInstanceOf[Call], callOutputsAsMap.map {
        case (k, v) => k -> v.getOrElse {
          throw new WdlExpressionException(s"Could not retrieve output $k value for call ${outputKey.scope.unqualifiedName}")
        }
      }))
    }
  }

  private def fetchCallInputEntries(callKey: ExecutionStoreKey): Traversable[SymbolStoreEntry] = {
    symbolCache.getOrElse(SymbolCacheKey(callKey.scope.fullyQualifiedName, input = true), Seq.empty)
  }

  /**
   * Load whatever execution statuses and symbols are stored for this workflow, regardless of whether this is a
   * workflow being restarted, or started for the first time.
   */
  private def createCaches: Future[(ExecutionStore, SymbolCache)] = {

    def isInScatterBlock(c: Call) = c.ancestry.exists(_.isInstanceOf[Scatter])
    import OldStyleFinalCall.FinalCallString

    val futureExecutionCache = globalDataAccess.getExecutionStatuses(workflow.id) map { statuses =>
      statuses map { case (k, v) =>
        val key: ExecutionStoreKey = if (k.fqn.isFinalCall) {
          // Final calls are not part of the workflow namespace, handle these differently from other keys.
          k.fqn.storeKey(workflow)
        } else {
          (workflow.namespace.resolve(k.fqn), k.index, k.attempt) match {
            case (Some(c: Call), Some(i), a) => BackendCallKey(c, Option(i), a)
            case (Some(c: Call), None, _) if isInScatterBlock(c) => CollectorKey(c)
            case (Some(c: Call), None, a) => BackendCallKey(c, None, a)
            case (Some(s: Scatter), None, _) => ScatterKey(s, None)
            case _ => throw new UnsupportedOperationException(s"Execution entry invalid: $k -> $v")
          }
        }
        key -> v.executionStatus
      }
    }

    val futureSymbolCache = globalDataAccess.getAllSymbolStoreEntries(workflow.id) map { symbols =>
      symbols groupBy(symbol => SymbolCacheKey(symbol.scope, symbol.isInput))
    }

    for {
      executionCache <- futureExecutionCache
      symbolCache <- futureSymbolCache
    } yield (executionCache, symbolCache)
  }

  private def buildSymbolStoreEntries(descriptor: OldStyleWorkflowDescriptor, inputs: HostInputs): Traversable[SymbolStoreEntry] = {
    val inputSymbols = inputs map {
      case (name, value) => SymbolStoreEntry(name, value, descriptor.hash(value), input = true)
    }

    val callSymbols = for {
      call <- descriptor.namespace.workflow.calls
      (k, v) <- call.inputMappings
    } yield SymbolStoreEntry(s"${call.fullyQualifiedName}.$k", v, descriptor.hash(v), input = true)

    inputSymbols.toSet ++ callSymbols.toSet
  }

  /**
   * This is the lookup function used to evaluate scatter collection expressions.
   *
   * For example, scatter(x in foo.bar) would evaluate the collection "foo.bar"
   * and call this lookup function on "foo".
   *
   * This implementation takes a few shortcuts and tries to find a Call or
   * Declaration with the given name in the workflow.
   *
   * A more long-term approach would be to traverse the scope hierarchy to resolve a variable into
   * the closest definition in scope.
   */
  private def scatterCollectionLookupFunction(workflow: Workflow, key: ExecutionStoreKey)(identifier: String): WdlValue = {
    resolveIdentifierOrElse(identifier, lookupCall(key, workflow), lookupDeclaration(workflow)) {
      throw new WdlExpressionException(s"Could not resolve identifier '$identifier' as a call or declaration.")
    }
  }

  private def isWorkflowDone: Boolean = executionStore forall executionStore.isDone

  private def isWorkflowAborted: Boolean = executionStore.values forall { state => isTerminal(state) || state == ExecutionStatus.NotStarted }

  private def processRunnableScatter(scatterKey: ScatterKey): Try[ExecutionStartResult] = {

    def buildExecutionStartResult(collection: WdlValue): ExecutionStartResult = {
      collection match {
        case a: WdlArray =>
          val newEntries = scatterKey.populate(a.value.size)
          executionStore += scatterKey -> ExecutionStatus.Starting
          executionStore ++= newEntries.keys map { _ -> ExecutionStatus.NotStarted }

          val persistFuture = for {
            _ <- persistStatus(scatterKey, ExecutionStatus.Starting, None)
            _ <- globalDataAccess.insertCalls(workflow.id, newEntries.keys, backend)
            _ = self ! PersistencesCompleted(newEntries.keys, ExecutionStatus.NotStarted)
          } yield ()

          persistFuture onComplete {
            case Success(_) => self ! ScatterCompleted(scatterKey)
            case Failure(t) =>
              logger.error(s"Error persisting status / inserting calls for scatter key ${scatterKey.tag} at status Starting", t)
              self ! CheckForWorkflowComplete
          }
          val startEntries = (newEntries.keys map { StartEntry(_, ExecutionStatus.NotStarted) } toSet) + StartEntry(scatterKey, ExecutionStatus.Starting)
          ExecutionStartResult(startEntries, newEntries = newEntries.keys)
        case v: WdlValue => throw new Throwable("Scatter collection must evaluate to an array")
      }
    }

    val rootWorkflow = scatterKey.scope.rootWorkflow
    for {
      collection <- scatterKey.scope.collection.evaluate(scatterCollectionLookupFunction(rootWorkflow, scatterKey), engineFunctions)
    } yield buildExecutionStartResult(collection)
  }

  private def processRunnableCollector(collector: CollectorKey): Try[ExecutionStartResult] = {
    executionStore += collector -> ExecutionStatus.Starting
    persistStatus(collector, ExecutionStatus.Starting)
    val shards: Iterable[BackendCallKey] = findShardEntries(collector) collect { case (k: BackendCallKey, v) if v == ExecutionStatus.Done => k }

    generateCollectorOutput(collector, shards) match {
      case Failure(e) =>
        self ! CallFailedNonRetryable(collector, Seq.empty, None, e.getMessage)
      case Success(outputs) =>
        logger.info(s"Collection complete for Scattered Call ${collector.tag}.")
        self ! CallCompleted(collector, outputs, Seq.empty, 0, hash = None, resultsClonedFrom = None)
    }

    Success(ExecutionStartResult(Set(StartEntry(collector, ExecutionStatus.Starting))))
  }

  private def sendStartMessage(callKey: BackendCallKey, callInputs: Map[String, WdlValue]) = {
    val descriptor = OldStyleBackendCallJobDescriptor(workflow, callKey, callInputs)
    val log = WorkflowLogger("WorkflowActor", workflow, akkaLogger = Option(akkaLogger))

    def startCallMsg = InitialStartCall(callKey, OldStyleCallActor.StartBackendCall(Option(backend.stdoutStderr(descriptor))))

    def loadCachedBackendCallAndMessage(workflow: OldStyleWorkflowDescriptor, cachedExecution: Execution) = {
      workflow.namespace.resolve(cachedExecution.callFqn) match {
        case Some(c: Call) =>
          val jobDescriptor = OldStyleBackendCallJobDescriptor(workflow, BackendCallKey(c, cachedExecution.index.toIndex, cachedExecution.attempt), callInputs)
          log.info(s"Call Caching: Cache hit. Using UUID(${jobDescriptor.workflowDescriptor.id.shortString}):${jobDescriptor.key.tag} as results for UUID(${descriptor.workflowDescriptor.id.shortString}):${descriptor.key.tag}")
          self ! UseCachedCall(callKey, OldStyleCallActor.UseCachedCall(jobDescriptor, descriptor,
            backend.stdoutStderr(descriptor)))
        case _ =>
          log.error(s"Call Caching: error when resolving '${cachedExecution.callFqn}' in workflow with execution ID ${cachedExecution.workflowExecutionId}: falling back to normal execution")
          self ! startCallMsg
      }
    }

    /* Tries to use the cached Execution to send a UseCachedCall message.  If anything fails, send an InitialStartCall message */
    def loadCachedCallOrInitiateCall(cachedDescriptor: Try[OldStyleWorkflowDescriptor], cachedExecution: Execution) = cachedDescriptor match {
      case Success(w) => loadCachedBackendCallAndMessage(w, cachedExecution)
      case Failure(ex) =>
        log.error(s"Call Caching: error when loading workflow with execution ID ${cachedExecution.workflowExecutionId}: falling back to normal execution", ex)
        self ! startCallMsg
    }

    def startCachedCall(cachedExecution: Execution)(implicit actorSystem: ActorSystem) = {
      val wfDesc = globalDataAccess.getWorkflowExecutionAndAux(cachedExecution.workflowExecutionId) flatMap workflowDescriptorFromExecutionAndAux
      wfDesc onComplete { cachedDescriptor =>
        loadCachedCallOrInitiateCall(cachedDescriptor, cachedExecution)
      }
    }

    def checkCacheAndStartCall = {
      log.debug(s"Call caching 'readFromCache' is turned on. Checking cache before starting call")
      descriptor.hash map { hash =>
        globalDataAccess.getExecutionsWithResuableResultsByHash(hash.overallHash) onComplete {
          case Success(executions) if executions.nonEmpty => startCachedCall(executions.head)
          case Success(_) =>
            log.info(s"Call Caching: cache miss")
            self ! startCallMsg
          case Failure(ex) =>
            log.error(s"Call Caching: Failed to look up executions that matched hash '$hash'. Falling back to normal execution", ex)
            self ! startCallMsg
        }
      } recover { case e =>
        log.error(s"Failed to calculate hash for call '${descriptor.key.tag}'.", e)
        self ! CallFailedToInitialize(callKey, s"Failed to calculate hash for call '${descriptor.key.tag}': ${e.getMessage}")
      }
    }

    def startCall = {
      if (descriptor.workflowDescriptor.readFromCache) {
        checkCacheAndStartCall
      } else {
        log.debug(s"Call caching 'readFromCache' is turned off, starting call")
        self ! startCallMsg
      }
    }

    Try(descriptor.callRuntimeAttributes) map { attrs =>
      val databaseKey = ExecutionDatabaseKey(descriptor.key.scope.fullyQualifiedName, descriptor.key.index, descriptor.key.attempt)
      globalDataAccess.upsertRuntimeAttributes(workflow.id, databaseKey, attrs.attributes, context.system) onComplete {
        case Success(_) => startCall
        case Failure(f) =>
          logger.error("Could not persist runtime attributes", f)
          self ! CallFailedToInitialize(callKey, s"Could not persist runtime attributes: ${f.getMessage}")
      }
    } recover {
      case f =>
        logger.error(f.getMessage, f)
        self ! CallFailedToInitialize(callKey, f.getMessage)
    }
  }

  private def processRunnableCall(callKey: CallKey): Try[ExecutionStartResult] = {
    // In the `startRunnableCalls` context, record the call as Starting and initiate persistence.
    // The restart scenario assumes a restartable/resumable call is already in Running.
    executionStore += callKey -> ExecutionStatus.Starting
    persistStatus(callKey, ExecutionStatus.Starting)

    callKey match {
      case backendCallKey: BackendCallKey =>
        fetchLocallyQualifiedInputs(backendCallKey) match {
          case Success(callInputs) =>
            sendStartMessage(backendCallKey, callInputs)
            Success(ExecutionStartResult(Set(StartEntry(callKey, ExecutionStatus.Starting))))
          case Failure(t) =>
            logger.error(s"Failed to fetch locally qualified inputs for call ${callKey.tag}", t)
            Failure(t)
        }
      case finalCallKey: FinalCallKey =>
        self ! InitialStartCall(finalCallKey, OldStyleCallActor.StartFinalCall)
        Success(ExecutionStartResult(Set(StartEntry(callKey, ExecutionStatus.Starting))))
    }
  }

  private def symbolsAsTable: Future[Traversable[Seq[String]]] = {
    def buildColumnData(entry: SymbolStoreEntry): Seq[String] = {
      // If the value of this symbol store entry is defined, that value truncated to MarkdownMaxColumnChars,
      // otherwise the empty string.
      def columnValue(entry: SymbolStoreEntry): String = {
        entry.wdlValue map { value =>
          s"(${value.wdlType.toWdlString}) ${value.valueString}".take(MarkdownMaxColumnChars)
        } getOrElse ""
      }
      // The Markdown columns as a Seq[String].
      Seq(
        entry.key.scope,
        entry.key.name,
        entry.key.index.map(_.toString).getOrElse(""),
        if (entry.key.input) "INPUT" else "OUTPUT",
        entry.wdlType.toWdlString,
        columnValue(entry)
      )
    }

    for {
      symbols <- globalDataAccess.getAllSymbolStoreEntries(workflow.id)
      columnData = symbols map buildColumnData
    } yield columnData
  }

  private def symbolsMarkdownTable: Future[Option[String]] = {
    val header = Seq("SCOPE", "NAME", "INDEX", "I/O", "TYPE", "VALUE")
    symbolsAsTable map {
      case rows if rows.isEmpty => None
      case rows => Option(TerminalUtil.mdTable(rows.toSeq, header))
    }
  }

  private def executionsAsTable: Future[Iterable[Seq[String]]] = {
    def buildColumnData(statusEntry: (ExecutionDatabaseKey, CallStatus)): Seq[String] = {
      val (k, v) = statusEntry
      Seq(k.fqn.toString, k.index.getOrElse("").toString, v.executionStatus.toString, v.returnCode.getOrElse("").toString)
    }

    for {
      executionStatuses <- globalDataAccess.getExecutionStatuses(workflow.id)
      columnData = executionStatuses map buildColumnData
    } yield columnData
  }

  private def executionsMarkdownTable: Future[Option[String]] = {
    val header = Seq("SCOPE", "INDEX", "STATUS", "RETURN CODE")
    executionsAsTable map {
      case rows if rows.isEmpty => None
      case rows => Option(TerminalUtil.mdTable(rows.toSeq, header))
    }
  }


  private def shutDown(): Unit = {
    logger.debug(s"WorkflowActor is done, shutting down.")
    context.stop(self)
  }

  override def postStop() = {
    super.postStop()
    /*
    NOTE: Due to the actor lifecycle, this may erase any older log files in the event of an actor restart
    as `Actor.preRestart` calls `Actor.postStop`.
     */
    workflow.maybeDeleteWorkflowLog()
  }
}
