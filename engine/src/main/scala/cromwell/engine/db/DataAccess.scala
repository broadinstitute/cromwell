package cromwell.engine.db

import java.sql.Timestamp
import java.util.Date

import akka.actor.ActorSystem
import cromwell.backend.{ExecutionEventEntry, ExecutionHash, JobKey}
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.core.{CallOutput, CallOutputs, WorkflowId}
import cromwell.database.SqlConverters._
import cromwell.database.SqlDatabase
import cromwell.database.obj._
import cromwell.database.slick.SlickDatabase
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus._
import cromwell.engine._
import cromwell.engine.backend.{OldStyleBackend, OldStyleBackendCallJobDescriptor, _}
import cromwell.engine.db.DataAccess.{RetryBackoff, WorkflowExecutionAndAux}
import cromwell.engine.db.EngineConverters._
import cromwell.engine.finalcall.OldStyleFinalCall
import cromwell.engine.workflow.OldStyleWorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.{BackendCallKey, ExecutionStoreKey, _}
import cromwell.services._
import cromwell.webservice.{CallCachingParameters, QueryMetadata, WorkflowQueryParameters, WorkflowQueryResponse}
import org.apache.commons.lang3.StringUtils
import org.joda.time.DateTime
import org.slf4j.LoggerFactory
import wdl4s.types.WdlPrimitiveType
import wdl4s.values.WdlValue
import wdl4s.{CallInputs, _}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

object DataAccess {
  lazy val log = LoggerFactory.getLogger(classOf[DataAccess])

  val globalDataAccess: DataAccess = new SlickDatabase() with DataAccess

  case class WorkflowExecutionAndAux(execution: WorkflowExecution, aux: WorkflowExecutionAux)

  val FailureEventMaxMessageLength = 1024
  val RetryBackoff = SimpleExponentialBackoff(50 millis, 1 seconds, 1D)
}

trait DataAccess extends AutoCloseable {
  this: SqlDatabase =>

  private def withRetry[A](f: () => Future[A])(implicit actorSystem: ActorSystem): Future[A] = {
    Retry.withRetry(f, maxRetries = Option(10), backoff = RetryBackoff, isTransient = isTransient)
  }

  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  def createWorkflow(workflowDescriptor: OldStyleWorkflowDescriptor,
                     workflowInputs: Traversable[SymbolStoreEntry],
                     calls: Traversable[Scope],
                     backend: OldStyleBackend)(implicit ec: ExecutionContext): Future[Unit] = {
    val workflowExecution = new WorkflowExecution(
      workflowDescriptor.id.toString,
      workflowDescriptor.name,
      WorkflowSubmitted.toString,
      new Date().toTimestamp,
      None)
    val workflowExecutionAux = (workflowExecutionId: Int) => {
      new WorkflowExecutionAux(
        workflowExecutionId,
        workflowDescriptor.sourceFiles.wdlSource.toClob,
        workflowDescriptor.sourceFiles.inputsJson.toClob,
        workflowDescriptor.sourceFiles.workflowOptionsJson.toClob
      )
    }
    val workflowSymbols = (workflowExecutionId: Int) => {
      workflowInputs.toSeq map {
        _.toSymbol(workflowDescriptor.namespace.workflow.outputs)(workflowExecutionId)
      }
    }
    val scopeKeys: Seq[ExecutionStoreKey] = calls.toSeq collect {
      case call: Call => BackendCallKey(call, None, 1)
      case scatter: Scatter => ScatterKey(scatter, None)
      case finalCall: OldStyleFinalCall => FinalCallKey(finalCall)
    }
    val executions = (workflowExecutionId: Int) => {
      scopeKeys map { key =>
        new Execution(
          workflowExecutionId = workflowExecutionId,
          callFqn = key.scope.fullyQualifiedName,
          index = key.index.fromIndex,
          attempt = key.attempt,
          status = ExecutionStatus.NotStarted.toString,
          rc = None,
          startDt = None,
          endDt = None,
          backendType = backend.backendType.displayName,
          allowsResultReuse = true,
          dockerImageHash = None,
          resultsClonedFrom = None,
          overallHash = None)
      }
    }
    val executionInfos = (executionId: Int) => {
      backend.executionInfoKeys map {
        ExecutionInfo(executionId, _, None)
      }
    }
    createWorkflow(workflowExecution, workflowExecutionAux, workflowSymbols, executions, executionInfos)
  }

  def getWorkflowState(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Option[WorkflowState]] = {
    getWorkflowStateString(workflowId.toString) map {
      _ map WorkflowState.fromString
    }
  }

  def assertWorkflowExistsByState(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = {
    getWorkflowState(workflowId) flatMap {
      case Some(_) => Future.successful(Unit)
      case None => Future.failed(new WorkflowNotFoundException(s"Workflow '$workflowId' not found"))
    }
  }

  def getWorkflowExecutionAndAux(workflowId: WorkflowId)
                                (implicit ec: ExecutionContext): Future[WorkflowExecutionAndAux] = {
    getWorkflowExecutionAndAuxTuple(workflowId.toString) map {
      case (execution, aux) => WorkflowExecutionAndAux(execution, aux)
    }
  }

  def getWorkflowExecutionAndAux(workflowExecutionId: Int)
                                (implicit ec: ExecutionContext): Future[WorkflowExecutionAndAux] = {
    getWorkflowExecutionAndAuxTuple(workflowExecutionId) map {
      case (execution, aux) => WorkflowExecutionAndAux(execution, aux)
    }
  }

  def getWorkflowExecutionAndAuxByState(states: Traversable[WorkflowState])
                                       (implicit ec: ExecutionContext): Future[Traversable[WorkflowExecutionAndAux]] = {
    getWorkflowExecutionAndAuxTuples(states.map(_.toString)) map {
      _ map { case (execution, aux) => WorkflowExecutionAndAux(execution, aux) }
    }
  }

  def getExecutionInfos(workflowId: WorkflowId, call: Call, attempt: Int)
                       (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]] = {
    getExecutionInfos(workflowId.toString, call.fullyQualifiedName, attempt)
  }

  def getExecutionInfoByKey(workflowId: WorkflowId, call: Call, attempt: Int, key: String)
                           (implicit ec: ExecutionContext): Future[Option[Option[String]]] = {
    getExecutionInfoByKey(workflowId.toString, call.fullyQualifiedName, attempt, key)
  }

  def updateExecutionInfo(workflowId: WorkflowId, callKey: BackendCallKey, key: String, value: Option[String])
                         (implicit ec: ExecutionContext): Future[Unit] = {
    updateExecutionInfo(workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex, callKey.attempt,
      key, value)
  }

  /**
    * TODO: the interface for retrying SQL commands that might fail because
    * of transient reasons (e.g. upsertExecutionInfo and upsertRuntimeAttributes)
    * could be made better.  Also it'd be nice if it were transparently
    * turned on/off for all methods in dataAccess.
    *
    * https://github.com/broadinstitute/cromwell/issues/693
    */
  def upsertExecutionInfo(workflowId: WorkflowId,
                          callKey: JobKey,
                          keyValues: Map[String, Option[String]],
                          actorSystem: ActorSystem): Future[Unit] = {
    implicit val system = actorSystem
    implicit val ec = actorSystem.dispatcher
    withRetry(() => upsertExecutionInfo(workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex,
      callKey.attempt, keyValues))
  }

  def upsertRuntimeAttributes(workflowId: WorkflowId,
                              key: ExecutionDatabaseKey,
                              attributes: Map[String, WdlValue],
                              actorSystem: ActorSystem): Future[Unit] = {
    implicit val system = actorSystem
    implicit val ec = actorSystem.dispatcher
    val keyValues = attributes mapValues { _.valueString }
    withRetry(() => upsertRuntimeAttributes(workflowId.toString, key.fqn, key.index.fromIndex, key.attempt, keyValues))
  }

  def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState)
                         (implicit ec: ExecutionContext): Future[Unit] = {
    updateWorkflowState(workflowId.toString, workflowState.toString,
      if (workflowState.isTerminal) Option(new Date().toTimestamp) else None
    )
  }

  def getAllSymbolStoreEntries(workflowId: WorkflowId)
                              (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val symbols = getAllSymbols(workflowId.toString)
    symbols map { _ map { _.toSymbolStoreEntry } }
  }

  // TODO needed to support compatibility with current code, this seems like an inefficient way of getting
  // TODO workflow outputs.
  /** Returns all outputs for this workflowId */
  def getWorkflowOutputs(workflowId: WorkflowId)
                        (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val symbols = getWorkflowOutputSymbols(workflowId.toString)
    symbols map { _ map { _.toSymbolStoreEntry } }
  }

  def getAllOutputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val symbols = getAllSymbols(workflowId.toString, IoOutput)
    symbols map { _ map { _.toSymbolStoreEntry } }
  }

  def getAllInputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val symbols = getAllSymbols(workflowId.toString, IoInput)
    symbols map { _ map { _.toSymbolStoreEntry } }
  }

  /** Get all outputs for the scope of this call. */
  def getOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val symbols = getAllSymbols(workflowId.toString, IoOutput, key.fqn, key.index.fromIndex)
    symbols map { _ map { _.toSymbolStoreEntry } }
  }

  def getAllRuntimeAttributes(id: WorkflowId)
                             (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, Map[String, String]]] = {
    val results = getAllRuntimeAttributes(id.toString)
    results map {
      _ map {
        case (callFqn, index, attempt, key, value) =>
          ExecutionDatabaseKey(callFqn, index.toIndex, attempt) -> (key -> value)
      } groupBy {
        case (key, _) => key
      } mapValues {
        _.map({ case (_, nameValue) => nameValue })
      } mapValues {
        _.toMap
      }
    }
  }

  /** Get all inputs for the scope of this call. */
  def getInputs(id: WorkflowId, call: Call)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val symbols = getAllSymbols(id.toString, IoInput, call.fullyQualifiedName, ExecutionIndex.IndexNone)
    symbols map { _ map { _.toSymbolStoreEntry } }
  }

  private def wdlValueToDbValue(v: WdlValue): String = v.wdlType match {
    case p: WdlPrimitiveType => v.valueString
    case o => v.toWdlString
  }

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  def setOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey, callOutputs: CallOutputs,
                 workflowOutputFqns: Seq[ReportableSymbol])(implicit ec: ExecutionContext): Future[Unit] = {
    val reportableResultNames = workflowOutputFqns map { _.fullyQualifiedName }
    val outputSymbols = (workflowExecutionId: Int) => {
      callOutputs.toSeq map {
        case (symbolLocallyQualifiedName, CallOutput(wdlValue, hash)) =>
          val reportableSymbol = key.index.fromIndex == -1 &&
            reportableResultNames.contains(key.fqn + "." + symbolLocallyQualifiedName)
          val value = wdlValueToDbValue(wdlValue).toNonEmptyClob
          new Symbol(
            workflowExecutionId,
            key.fqn,
            symbolLocallyQualifiedName,
            key.index.fromIndex,
            IoOutput,
            reportableSymbol,
            wdlValue.wdlType.toWdlString,
            value,
            hash.map(_.value)
          )
      }
    }
    setOutputs(workflowId.toString, outputSymbols)
  }

  /** Updates the existing input symbols to replace expressions with real values **/
  def updateCallInputs(workflowId: WorkflowId, key: BackendCallKey, callInputs: CallInputs)
                      (implicit ec: ExecutionContext): Future[Int] = {

    val mappedInputs = callInputs.toSeq map {
      case (inputName, wdlValue) =>
        (inputName, wdlValue.wdlType.toWdlString, Option(wdlValueToDbValue(wdlValue).toClob))
    }
    updateCallInputs(workflowId.toString, key.scope.fullyQualifiedName, key.index.fromIndex, mappedInputs).map(_.sum)
  }

  def setExecutionEvents(workflowId: WorkflowId, callFqn: String, shardIndex: Option[Int], attempt: Int,
                         events: Seq[ExecutionEventEntry])(implicit ec: ExecutionContext): Future[Unit] = {
    val executionEvents = (executionId: Int) => {
      events map { executionEventEntry =>
        new ExecutionEvent(
          executionId,
          executionEventEntry.description,
          new Timestamp(executionEventEntry.startTime.getMillis),
          new Timestamp(executionEventEntry.endTime.getMillis))
      }
    }
    shardIndex match {
      case Some(idx) =>
        setExecutionEvents(workflowId.toString, callFqn, idx, attempt, executionEvents)
      case None =>
        setExecutionEvents(workflowId.toString, callFqn, attempt, executionEvents)
    }
  }

  /** Gets a mapping from call FQN to an execution event entry list */
  def getAllExecutionEvents(workflowId: WorkflowId)(implicit ec: ExecutionContext):
  Future[Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]]] = {
    // The database query gives us a Seq[(CallFqn, ExecutionEvent)]. We want a Map[CallFqn -> ExecutionEventEntry].
    // So let's do some functional programming!
    val results = getAllExecutionEvents(workflowId.toString)
    results map toExecutionEvents
  }

  private def toExecutionEvents(events: Traversable[(String, Int, Int, String, Timestamp, Timestamp)])
                               (implicit ec: ExecutionContext): Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]] = {
    val mapped: Traversable[(ExecutionDatabaseKey, ExecutionEventEntry)] = events map {
      case (fqn, index, attempt, description, startTime, endTime) =>
        ExecutionDatabaseKey(fqn, index.toIndex, attempt) ->
          ExecutionEventEntry(description, new DateTime(startTime), new DateTime(endTime))
    }
    val grouped: Map[ExecutionDatabaseKey, Seq[(ExecutionDatabaseKey, ExecutionEventEntry)]] = mapped.toSeq groupBy {
      case (key, _) => key
    }
    grouped mapValues {
      _ map { case (_, values) => values }
    }
  }

  /** Add a new failure event for a call into the database. */
  def addCallFailureEvent(workflowId: WorkflowId, executionKey: ExecutionDatabaseKey,
                          failure: FailureEventEntry)(implicit ec: ExecutionContext): Future[Unit] = {
    val failureEvents = (workflowExecutionId: Int, executionId: Int) => {
      new FailureEvent(
        workflowExecutionId,
        Option(executionId),
        StringUtils.abbreviate(failure.failure, DataAccess.FailureEventMaxMessageLength),
        new Timestamp(failure.timestamp.getMillis))
    }
    executionKey match {
      case ExecutionDatabaseKey(callFqn, Some(index), attempt) =>
        addCallFailureEvent(workflowId.toString, callFqn, index, attempt, failureEvents)
      case ExecutionDatabaseKey(callFqn, None, attempt) =>
        addCallFailureEvent(workflowId.toString, callFqn, attempt, failureEvents)
    }
  }

  /** Add a new failure event for a workflow into the database. */
  def addWorkflowFailureEvent(workflowId: WorkflowId, failure: FailureEventEntry)
                             (implicit ec: ExecutionContext): Future[Unit] = {

    val failureEvents = (workflowExecutionId: Int) => {
      new FailureEvent(
        workflowExecutionId,
        None,
        StringUtils.abbreviate(failure.failure, DataAccess.FailureEventMaxMessageLength),
        new Timestamp(failure.timestamp.getMillis))
    }
    addWorkflowFailureEvent(workflowId.toString, failureEvents)
  }

  def addMetadataEvent(metadataEvent: MetadataEvent)(implicit ec: ExecutionContext): Future[Unit] = {
    val key = metadataEvent.key
    val workflowUuid = key.workflowId.id.toString
    val timestamp = metadataEvent.timestamp
    key.jobKey match {
      case None => addMetadataEvent(workflowUuid, key.key, metadataEvent.value.value, timestamp)
      case Some(jobKey) => addMetadataEvent(workflowUuid, key.key, jobKey.callFqn, jobKey.index, jobKey.attempt, metadataEvent.value.value, timestamp)
    }
  }

  def queryMetadataEvents(query: MetadataQuery)(implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = query.workflowId.id.toString
    val futureMetadata: Future[Seq[Metadatum]] = query match {
      case MetadataQuery(_, None, None) => queryMetadataEvents(uuid)
      case MetadataQuery(_, None, Some(key)) => queryMetadataEvents(uuid, key)
      case MetadataQuery(_, Some(jobKey), None) => queryMetadataEvents(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, Some(jobKey), Some(key)) => queryMetadataEvents(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt)
    }
    futureMetadata map { metadata =>
      metadata map { m =>
        // If callFqn is non-null then attempt will also be non-null and there is a MetadataJobKey.
        val metadataJobKey: Option[MetadataJobKey] = for {
          callFqn <- m.callFqn
          attempt <- m.attempt
        } yield new MetadataJobKey(callFqn, m.index, attempt)

        val key = MetadataKey(query.workflowId, metadataJobKey, m.key)
        val value = MetadataValue(m.value.orNull)
        MetadataEvent(key, value, m.timestamp)
      }
    }
  }

  /** Retrieve all recorded failures for a Workflow */
  def getFailureEvents(workflowId: WorkflowId)
                      (implicit ec: ExecutionContext): Future[Seq[QualifiedFailureEventEntry]] = {
    val results = getFailureEvents(workflowId.toString)
    results map {
      _.toSeq map {
        case (workflowUuid, message, timestamp, Some(fqn), Some(index), Some(attempt)) =>
          val key = ExecutionDatabaseKey(fqn, index.toIndex, attempt)
          QualifiedFailureEventEntry(workflowUuid, Option(key), message, new DateTime(timestamp))
        case (workflowUuid, message, timestamp, None, None, None) =>
          QualifiedFailureEventEntry(workflowUuid, None, message, new DateTime(timestamp))
        case unexpected =>
          throw new RuntimeException(s"Unexpected failure event: $unexpected")
      }
    }
  }


  /** Set the status of one or several calls to starting and update the start date. */
  def setStartingStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey])
                       (implicit ec: ExecutionContext): Future[Unit] = {
    if (scopeKeys.isEmpty) {
      Future.successful(())
    } else {
      val mappedKeys = scopeKeys map { scopeKey =>
        (scopeKey.fqn, scopeKey.index.fromIndex, scopeKey.attempt)
      }
      setStartingStatus(workflowId.toString, ExecutionStatus.Starting.toString, Option(new Date().toTimestamp),
        mappedKeys)
    }
  }

  /** Simply set the status of one of several calls. Status cannot be Starting or a Terminal status. */
  def updateStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey], status: ExecutionStatus)
                  (implicit ec: ExecutionContext): Future[Unit] = {
    if (scopeKeys.isEmpty) {
      Future.successful(())
    } else {
      require(!status.isTerminal && !(status == ExecutionStatus.Starting))
      val mappedKeys = scopeKeys map { scopeKey =>
        (scopeKey.fqn, scopeKey.index.fromIndex, scopeKey.attempt)
      }
      updateStatus(workflowId.toString, status.toString, mappedKeys)
    }
  }

  /** Set the status of a Call to a terminal status, and update associated information (return code, hash, cache). */
  def setTerminalStatus(workflowId: WorkflowId, scopeKeys: ExecutionDatabaseKey, status: ExecutionStatus,
                        scriptReturnCode: Option[Int], hash: Option[ExecutionHash],
                        resultsClonedFrom: Option[OldStyleBackendCallJobDescriptor])
                       (implicit ec: ExecutionContext): Future[Unit] = {
    require(status.isTerminal)

    val workflowUuid = workflowId.toString
    val callFqn = scopeKeys.fqn
    val index = scopeKeys.index.fromIndex
    val attempt = scopeKeys.attempt
    val statusString = status.toString
    val endDt = Option(new Date().toTimestamp)
    val overallHash = hash map { _.overallHash }
    val dockerHash = hash flatMap { _.dockerHash }
    resultsClonedFrom match {
      case Some(jobDescriptor) =>
        val workflowUuidClone = jobDescriptor.workflowDescriptor.id.toString
        val callFqnClone = jobDescriptor.key.scope.fullyQualifiedName
        val indexClone = jobDescriptor.key.index.fromIndex
        val attemptClone = jobDescriptor.key.attempt
        setTerminalWithClone(workflowUuid, callFqn, index, attempt,
          statusString, endDt, scriptReturnCode, overallHash, dockerHash,
          workflowUuidClone, callFqnClone, indexClone, attemptClone)
      case None =>
        setTerminalWithoutClone(workflowUuid, callFqn, index, attempt,
          statusString, endDt, scriptReturnCode, overallHash, dockerHash)
    }
  }

  def getExecutionStatuses(workflowId: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]] = {
    val results = getExecutionStatuses(workflowId.toString)
    results map {
      _ map {
        case (callFqn, index, attempt, status, rc, hash, dockerHash) =>
          ExecutionDatabaseKey(callFqn, index.toIndex, attempt) ->
            CallStatus(status.toExecutionStatus, rc, hash map {
              ExecutionHash(_, dockerHash)
            }, None)
      } groupBy {
        case (key, _) => key
      } mapValues {
        _.map({ case (_, callStatus) => callStatus })
      } mapValues {
        _.head
      }
    }
  }

  /** Return all execution entries for the FQN, including collector and shards if any */
  def getExecutionStatuses(workflowId: WorkflowId, fqn: FullyQualifiedName)
                          (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]] = {
    val results = getExecutionStatuses(workflowId.toString, fqn)
    results map {
      _ map {
        case (callFqn, index, attempt, status, rc, hash, dockerHash) =>
          ExecutionDatabaseKey(callFqn, index.toIndex, attempt) ->
            CallStatus(status.toExecutionStatus, rc, hash map {
              ExecutionHash(_, dockerHash)
            }, None)
      } groupBy {
        case (key, _) => key
      } mapValues {
        _.map({ case (_, callStatus) => callStatus })
      } mapValues {
        _.head
      }
    }
  }

  def getExecutionStatus(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                        (implicit ec: ExecutionContext): Future[Option[CallStatus]] = {
    val results = getExecutionStatus(workflowId.toString, key.fqn, key.index.fromIndex, key.attempt)
    results map {
      _ map {
        case (status, rc, hash, dockerHash) =>
          CallStatus(status.toExecutionStatus, rc, hash map {
            ExecutionHash(_, dockerHash)
          }, None)
      }
    }
  }

  def insertCalls(workflowId: WorkflowId, keys: Traversable[ExecutionStoreKey], backend: OldStyleBackend)
                 (implicit ec: ExecutionContext): Future[Unit] = {
    val executions = (workflowExecutionId: Int) => {
      keys.toSeq map { key =>
        new Execution(
          workflowExecutionId = workflowExecutionId,
          callFqn = key.scope.fullyQualifiedName,
          index = key.index.fromIndex,
          attempt = key.attempt,
          status = ExecutionStatus.NotStarted.toString,
          rc = None,
          startDt = None,
          endDt = None,
          backendType = backend.backendType.displayName,
          allowsResultReuse = true,
          dockerImageHash = None,
          resultsClonedFrom = None,
          overallHash = None)
      }
    }
    val executionInfos = (executionId: Int) => {
      backend.executionInfoKeys map {
        ExecutionInfo(executionId, _, None)
      }
    }
    insertCalls(workflowId.toString, executions, executionInfos)
  }

  def getExecutionsAsExecutionInfos(id: WorkflowId)
                                   (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]] = {
    getExecutions(id) map { _ map { ExecutionInfosByExecution(_, Seq.empty) } }
  }

  def getExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[Execution]] = {
    getExecutions(id.toString)
  }

  override def getExecutionsWithResuableResultsByHash(hash: String)
                                                     (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  def updateWorkflowOptions(workflowId: WorkflowId, workflowOptionsJson: String)
                           (implicit ec: ExecutionContext): Future[Unit] = {
    updateWorkflowOptions(workflowId.toString, workflowOptionsJson.toClob)
  }

  def queryWorkflows(queryParameters: WorkflowQueryParameters)
                    (implicit ec: ExecutionContext): Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
    val workflowExecutions = queryWorkflowExecutions(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toDate.toTimestamp), queryParameters.endDate.map(_.toDate.toTimestamp),
      queryParameters.page, queryParameters.pageSize)

    val workflowCount = countWorkflowExecutions(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toDate.toTimestamp), queryParameters.endDate.map(_.toDate.toTimestamp))

    workflowCount flatMap { count =>
      workflowExecutions map { workflows =>
        (WorkflowQueryResponse(workflows.toSeq map { workflow =>
          WorkflowQueryResult(
            id = workflow.workflowExecutionUuid,
            name = workflow.name,
            status = workflow.status,
            start = new DateTime(workflow.startDt),
            end = workflow.endDt map {
              new DateTime(_)
            })
        }),
        //only return metadata if page is defined
        queryParameters.page map { _ => QueryMetadata(queryParameters.page, queryParameters.pageSize, Option(count)) })
      }
    }
  }

  def updateCallCaching(cachingParameters: CallCachingParameters)(implicit ec: ExecutionContext): Future[Int] = {
    // Figure out which of the three possible queries to use based on whether a call has been specified and
    // if so whether an index has been specified.
    (cachingParameters.callKey, cachingParameters.callKey flatMap {
      _.index
    }) match {
      case (Some(key), Some(idx)) =>
        updateCallCaching(cachingParameters.workflowId.toString, cachingParameters.allow, key.fqn, idx, key.attempt)
      case (Some(key), None) =>
        updateCallCaching(cachingParameters.workflowId.toString, cachingParameters.allow, key.fqn, key.attempt)
      case _ => updateCallCaching(cachingParameters.workflowId.toString, cachingParameters.allow)
    }
  }

  def infosByExecution(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]] = {
    infosByExecution(id.toString) map (_.toSeq) map ExecutionInfosByExecution.fromRawTuples
  }

  def infosByExecution(id: WorkflowId, fqn: FullyQualifiedName)
                      (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]] = {
    infosByExecution(id.toString, fqn) map (_.toSeq) map db.ExecutionInfosByExecution.fromRawTuples
  }

  /** Used by restart workflow code to get all executions that were in-flight when the engine is restarted. */
  def runningExecutionsAndExecutionInfos(id: WorkflowId)(implicit ec: ExecutionContext):
  Future[Traversable[ExecutionInfosByExecution]] = {
    val statuses = Set(ExecutionStatus.Starting, ExecutionStatus.Running) map (_.toString)
    val results = runningExecutionsAndExecutionInfos(id.toString, statuses)
    results map (_.toSeq) map db.ExecutionInfosByExecution.fromRawTuples
  }

  def callCacheDataByExecution(id: WorkflowId)
                              (implicit ec: ExecutionContext): Future[Traversable[ExecutionWithCacheData]] = {
    val tuple3s = callCacheDataByExecution(id.toString)
    tuple3s map {
      _.map {
        case (execution, cacheHitWorkflowUuid, cacheHitCallFqn) =>
          val cacheHit = for {
            workflowUuid <- cacheHitWorkflowUuid
            callFqn <- cacheHitCallFqn
          } yield CallCacheHit(workflowUuid, callFqn)

          ExecutionWithCacheData(execution, cacheHit)
      }
    }
  }
}
