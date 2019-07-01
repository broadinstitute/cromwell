package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import cats.syntax.apply._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.core.{HogGroup, WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreAbortResponse.WorkflowStoreAbortResponse
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{NotInOnHoldStateException, WorkflowStoreAbortResponse, WorkflowStoreState, WorkflowSubmissionResponse}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.collection._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

object SqlWorkflowStore {
  case class WorkflowSubmissionResponse(state: WorkflowStoreState, id: WorkflowId)

  case class NotInOnHoldStateException(workflowId: WorkflowId) extends
    Exception(
      s"Couldn't change status of workflow $workflowId to " +
        "'Submitted' because the workflow is not in 'On Hold' state"
    )

  object WorkflowStoreState extends Enumeration {
    type WorkflowStoreState = Value
    val Submitted = Value("Submitted")
    val Running = Value("Running")
    val Aborting = Value("Aborting")
    val OnHold = Value("On Hold")
  }

  object WorkflowStoreAbortResponse extends Enumeration {
    type WorkflowStoreAbortResponse = Value
    val NotFound = Value("NotFound")
    val AbortedOnHoldOrSubmitted = Value("AbortedOnHoldOrSubmitted")
    val AbortRequested = Value("AbortRequested")
  }
}

case class SqlWorkflowStore(sqlDatabase: WorkflowStoreSqlDatabase) extends WorkflowStore {
  /** This is currently hardcoded to success but used to do stuff, left in place for now as a useful
    *  startup initialization hook. */
  override def initialize(implicit ec: ExecutionContext): Future[Unit] = Future.successful(())

  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStoreAbortResponse] = {
    sqlDatabase.deleteOrUpdateWorkflowToState(
      workflowExecutionUuid = id.toString,
      workflowStateToDelete1 = WorkflowStoreState.OnHold.toString,
      workflowStateToDelete2 = WorkflowStoreState.Submitted.toString,
      workflowStateForUpdate = WorkflowStoreState.Aborting.toString
    ) map {
      case Some(true) =>
        WorkflowStoreAbortResponse.AbortedOnHoldOrSubmitted
      case Some(false) =>
        WorkflowStoreAbortResponse.AbortRequested
      case None =>
        WorkflowStoreAbortResponse.NotFound
    }
  }

  override def findWorkflows(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[WorkflowId]] = {
    sqlDatabase.findWorkflows(cromwellId) map { _ map WorkflowId.fromString }
  }

  override def findWorkflowsWithAbortRequested(cromwellId: String)(implicit ec: ExecutionContext): Future[Iterable[WorkflowId]] = {
    sqlDatabase.findWorkflowsWithAbortRequested(cromwellId) map { _ map WorkflowId.fromString }
  }

  override def abortAllRunning()(implicit ec: ExecutionContext): Future[Unit] = {
    sqlDatabase.setStateToState(WorkflowStoreState.Running.toString, WorkflowStoreState.Aborting.toString)
  }

  override def stats(implicit ec: ExecutionContext): Future[Map[WorkflowStoreState, Int]] = {
    sqlDatabase.workflowStateCounts.map {
      _ map {
        case (key, value) => WorkflowStoreState.withName(key) -> value
      }
    }
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    import cats.instances.list._
    import cats.syntax.traverse._
    import common.validation.Validation._
    sqlDatabase.fetchWorkflowsInState(
      n,
      cromwellId,
      heartbeatTtl.ago,
      OffsetDateTime.now.toSystemTimestamp,
      WorkflowStoreState.Submitted.toString,
      WorkflowStoreState.Running.toString,
      WorkflowStoreState.OnHold.toString
    ) map {
      // .get on purpose here to fail the future if something went wrong
      _.toList.traverse(fromWorkflowStoreEntry).toTry.get
    }
  }

  override def writeWorkflowHeartbeats(workflowIds: Set[(WorkflowId, OffsetDateTime)],
                                       heartbeatDateTime: OffsetDateTime)
                                      (implicit ec: ExecutionContext): Future[Int] = {
    val sortedWorkflowIds = workflowIds.toList sortBy(_._2) map (_._1.toString)
    sqlDatabase.writeWorkflowHeartbeats(sortedWorkflowIds, heartbeatDateTime.toSystemTimestamp)
  }

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowSubmissionResponse]] = {

    val asStoreEntries = sources map toWorkflowStoreEntry
    val returnValue = asStoreEntries map { workflowStore =>
      WorkflowSubmissionResponse(
        WorkflowStoreState.withName(workflowStore.workflowState),
        WorkflowId.fromString(workflowStore.workflowExecutionUuid)
      )
    }

    // The results from the Future aren't useful, so on completion map it into the precalculated return value instead. Magic!
    sqlDatabase.addWorkflowStoreEntries(asStoreEntries.toList) map { _ => returnValue }
  }

  override def switchOnHoldToSubmitted(id: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = {
    for {
      updated <- sqlDatabase.updateWorkflowState(
        id.toString,
        WorkflowStoreState.OnHold.toString,
        WorkflowStoreState.Submitted.toString
      )
      _ <- if (updated == 0) Future.failed(NotInOnHoldStateException(id)) else Future.successful(())
    } yield ()
  }


  private def fromWorkflowStoreEntry(workflowStoreEntry: WorkflowStoreEntry): ErrorOr[WorkflowToStart] = {

    val workflowOptionsValidation: ErrorOr[WorkflowOptions] = WorkflowOptions.fromJsonString(workflowStoreEntry.workflowOptions.toRawString).toErrorOr
    val startableStateValidation = workflowStoreStateToStartableState(workflowStoreEntry)

    (startableStateValidation, workflowOptionsValidation) mapN { (startableState, workflowOptions) =>

      val sources = WorkflowSourceFilesCollection(
        workflowSource = workflowStoreEntry.workflowDefinition.toRawStringOption,
        workflowUrl = workflowStoreEntry.workflowUrl,
        workflowRoot = workflowStoreEntry.workflowRoot,
        workflowType = workflowStoreEntry.workflowType,
        workflowTypeVersion = workflowStoreEntry.workflowTypeVersion,
        inputsJson = workflowStoreEntry.workflowInputs.toRawString,
        workflowOptions = workflowOptions,
        labelsJson = workflowStoreEntry.customLabels.toRawString,
        importsFile = workflowStoreEntry.importsZip.toBytesOption,
        warnings = Vector.empty,
        workflowOnHold = false
      )

      val id = WorkflowId.fromString(workflowStoreEntry.workflowExecutionUuid)
      val hogGroup: HogGroup = workflowStoreEntry.hogGroup.map(HogGroup(_)).getOrElse(HogGroup.decide(workflowOptions, id))

      WorkflowToStart(
        id = id,
        submissionTime = workflowStoreEntry.submissionTime.toSystemOffsetDateTime,
        sources = sources,
        state = startableState,
        hogGroup = hogGroup
      )
    }
  }

  private def workflowSubmissionState(workflowSourceFiles: WorkflowSourceFilesCollection) = {
    if (workflowSourceFiles.workflowOnHold)
      WorkflowStoreState.OnHold
    else
      WorkflowStoreState.Submitted
  }

  private def toWorkflowStoreEntry(workflowSourceFiles: WorkflowSourceFilesCollection): WorkflowStoreEntry = {
    import eu.timepit.refined._
    val nonEmptyJsonString: String Refined NonEmpty  = refineMV[NonEmpty]("{}")

    val actualWorkflowState = workflowSubmissionState(workflowSourceFiles)

    val workflowId = WorkflowId.randomId()
    val hogGroup = HogGroup.decide(workflowSourceFiles.workflowOptions, workflowId)

    WorkflowStoreEntry(
      workflowExecutionUuid = workflowId.toString,
      workflowDefinition = workflowSourceFiles.workflowSource.toClobOption,
      workflowUrl = workflowSourceFiles.workflowUrl,
      workflowRoot = workflowSourceFiles.workflowRoot,
      workflowType = workflowSourceFiles.workflowType,
      workflowTypeVersion = workflowSourceFiles.workflowTypeVersion,
      workflowInputs = workflowSourceFiles.inputsJson.toClobOption,
      workflowOptions = workflowSourceFiles.workflowOptions.asPrettyJson.toClobOption,
      customLabels = workflowSourceFiles.labelsJson.toClob(default = nonEmptyJsonString),
      workflowState = actualWorkflowState.toString,
      cromwellId = None,
      heartbeatTimestamp = None,
      submissionTime = OffsetDateTime.now.toSystemTimestamp,
      importsZip = workflowSourceFiles.importsZipFileOption.toBlobOption,
      hogGroup = Option(hogGroup.value)
    )
  }

  private def workflowStoreStateToStartableState(workflowStoreEntry: WorkflowStoreEntry): ErrorOr[StartableState] = {
    val workflowState = WorkflowStoreState.withName(workflowStoreEntry.workflowState)
    import cats.syntax.validated._
    // A workflow is startable if it's in Submitted, Running or Aborting state.
    workflowState match {
      case WorkflowStoreState.Submitted => Submitted.validNel
      case WorkflowStoreState.Running => RestartableRunning.validNel
      case WorkflowStoreState.Aborting => RestartableAborting.validNel
      case _ =>
        "Workflow %s in state %s cannot be started and should not have been fetched."
          .format(workflowStoreEntry.workflowExecutionUuid, workflowState)
          .invalidNel
    }
  }
}
