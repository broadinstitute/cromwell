package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import common.validation.ErrorOr.ErrorOr
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowSubmissionResponse
import eu.timepit.refined.api.Refined
import eu.timepit.refined.collection._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

object SqlWorkflowStore {
  case class WorkflowSubmissionResponse(state: WorkflowStoreState, id: WorkflowId)
}

case class SqlWorkflowStore(sqlDatabase: WorkflowStoreSqlDatabase) extends WorkflowStore {
  lazy val cromwellId = ConfigFactory.load().as[Option[String]]("system.cromwell_id")

  /** This is currently hardcoded to success but used to do stuff, left in place for now as a useful
    *  startup initialization hook. */
  override def initialize(implicit ec: ExecutionContext): Future[Unit] = Future.successful(())

  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    sqlDatabase.setToAborting(id.toString)
  }
  
  override def abortAllRunning()(implicit ec: ExecutionContext): Future[Unit] = {
    sqlDatabase.setAllRunningToAborting()
  }

  override def stats(implicit ec: ExecutionContext): Future[Map[WorkflowStoreState, Int]] = sqlDatabase.stats

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    sqlDatabase.removeWorkflowStoreEntry(id.toString).map(_ > 0) // i.e. did anything get deleted
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchStartableWorkflows(n: Int, cromwellId: String, heartbeatTtl: FiniteDuration)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    import cats.instances.list._
    import cats.syntax.traverse._
    import common.validation.Validation._
    sqlDatabase.fetchStartableWorkflows(n, cromwellId, heartbeatTtl) map {
      // .get on purpose here to fail the future if something went wrong
      _.toList.traverse(fromWorkflowStoreEntry).toTry.get
    }
  }

  override def writeWorkflowHeartbeats(workflowIds: Set[WorkflowId])(implicit ec: ExecutionContext): Future[Int] = {
    sqlDatabase.writeWorkflowHeartbeats(workflowIds map { _.toString })
  }

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowSubmissionResponse]] = {

    val asStoreEntries = sources map toWorkflowStoreEntry
    val returnValue = asStoreEntries map { workflowStore =>
      WorkflowSubmissionResponse(
        workflowStore.workflowState, WorkflowId.fromString(workflowStore.workflowExecutionUuid)
      )
    }

    // The results from the Future aren't useful, so on completion map it into the precalculated return value instead. Magic!
    sqlDatabase.addWorkflowStoreEntries(asStoreEntries.toList) map { _ => returnValue }
  }

  override def switchOnHoldToSubmitted(id: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = {
    sqlDatabase.setOnHoldToSubmitted(id.toString)
  }


  private def fromWorkflowStoreEntry(workflowStoreEntry: WorkflowStoreEntry): ErrorOr[WorkflowToStart] = {
    val sources = WorkflowSourceFilesCollection(
      workflowSource = workflowStoreEntry.workflowDefinition.toRawString,
      workflowRoot = workflowStoreEntry.workflowRoot,
      workflowType = workflowStoreEntry.workflowType,
      workflowTypeVersion = workflowStoreEntry.workflowTypeVersion,
      inputsJson = workflowStoreEntry.workflowInputs.toRawString,
      workflowOptionsJson = workflowStoreEntry.workflowOptions.toRawString,
      labelsJson = workflowStoreEntry.customLabels.toRawString,
      importsFile = workflowStoreEntry.importsZip.toBytesOption,
      warnings = Vector.empty,
      workflowOnHold = false
    )

    workflowStoreStateToStartableState(workflowStoreEntry) map { startableState =>
      WorkflowToStart(
        WorkflowId.fromString(workflowStoreEntry.workflowExecutionUuid),
        sources,
        startableState)
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

    WorkflowStoreEntry(
      workflowExecutionUuid = WorkflowId.randomId().toString,
      workflowDefinition = workflowSourceFiles.workflowSource.toClobOption,
      workflowRoot = workflowSourceFiles.workflowRoot,
      workflowType = workflowSourceFiles.workflowType,
      workflowTypeVersion = workflowSourceFiles.workflowTypeVersion,
      workflowInputs = workflowSourceFiles.inputsJson.toClobOption,
      workflowOptions = workflowSourceFiles.workflowOptionsJson.toClobOption,
      customLabels = workflowSourceFiles.labelsJson.toClob(default = nonEmptyJsonString),
      workflowState = actualWorkflowState,
      cromwellId = None,
      heartbeatTimestamp = None,
      submissionTime = OffsetDateTime.now.toSystemTimestamp,
      importsZip = workflowSourceFiles.importsZipFileOption.toBlobOption
    )
  }

  private def workflowStoreStateToStartableState(workflowStoreEntry: WorkflowStoreEntry): ErrorOr[StartableState] = {
    val workflowState = workflowStoreEntry.workflowState
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
