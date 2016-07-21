package cromwell.engine.workflow.workflowstore

import java.sql.Timestamp
import java.time.LocalDateTime

import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.database.WorkflowStoreSqlDatabase
import cromwell.database.obj.{WorkflowStoreEntry, WorkflowStoreEntryState}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

import scala.concurrent.{ExecutionContext, Future}
import scalaz.NonEmptyList

case class SqlWorkflowStore(sqlDatabase: WorkflowStoreSqlDatabase) extends WorkflowStore {
  override def initialize(implicit ec: ExecutionContext): Future[Unit] = sqlDatabase.initialize

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = sqlDatabase.remove(id.toString)

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchRunnableWorkflows(n: Int, state: StartableState)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    val dbState = state match {
      case WorkflowStoreState.Restartable => WorkflowStoreEntryState.Restartable
      case WorkflowStoreState.Submitted => WorkflowStoreEntryState.Submitted
    }
    sqlDatabase.fetchRunnableWorkflows(n, dbState) map { _ map fromWorkflowStoreEntry }
  }

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFiles])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowId]] = {

    val asStoreEntries = sources map toWorkflowStoreEntry
    val returnValue = asStoreEntries map { x => WorkflowId.fromString(x.workflowUuid) }

    // The results from the Future aren't useful, so on completion map it into the precalculated return value instead. Magic!
    sqlDatabase.add(sources.map(toWorkflowStoreEntry).list) map { _ => returnValue }
  }

  private def fromWorkflowStoreEntry(workflowStoreEntry: WorkflowStoreEntry): WorkflowToStart = {
    val sources = WorkflowSourceFiles(
      workflowStoreEntry.workflowSource,
      workflowStoreEntry.workflowInputs.getOrElse("{}"),
      workflowStoreEntry.workflowOptions.getOrElse("{}"))
    WorkflowToStart(
      WorkflowId.fromString(workflowStoreEntry.workflowUuid),
      sources,
      fromDbStateStringToStartableState(workflowStoreEntry.state))
  }

  private def toWorkflowStoreEntry(workflowSourceFiles: WorkflowSourceFiles): WorkflowStoreEntry = {
    WorkflowStoreEntry(
      WorkflowId.randomId().toString,
      workflowSourceFiles.wdlSource,
      Option(workflowSourceFiles.inputsJson),
      Option(workflowSourceFiles.workflowOptionsJson),
      WorkflowStoreEntryState.Submitted.toString,
      timestamp = Timestamp.valueOf(LocalDateTime.now()),
      workflowStoreId = None
    )
  }

  private def fromDbStateStringToStartableState(dbStateName: String): StartableState = {
    if (dbStateName.equalsIgnoreCase(WorkflowStoreEntryState.Restartable.toString)) {
      WorkflowStoreState.Restartable
    }
    else {
      WorkflowStoreState.Submitted
    }
  }
}
