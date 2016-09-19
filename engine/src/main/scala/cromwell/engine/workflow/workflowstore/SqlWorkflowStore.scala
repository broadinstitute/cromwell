package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

import scala.concurrent.{ExecutionContext, Future}
import scalaz.NonEmptyList

case class SqlWorkflowStore(sqlDatabase: WorkflowStoreSqlDatabase) extends WorkflowStore {
  override def initialize(implicit ec: ExecutionContext): Future[Unit] = {
    sqlDatabase.updateWorkflowState(
      WorkflowStoreState.Running.toString,
      WorkflowStoreState.Restartable.toString)
  }

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    sqlDatabase.removeWorkflowStoreEntry(id.toString).map(_ > 0) // i.e. did anything get deleted
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchRunnableWorkflows(n: Int, state: StartableState)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    sqlDatabase.queryWorkflowStoreEntries(n, state.toString, WorkflowStoreState.Running.toString) map {
      _.toList map fromWorkflowStoreEntry
    }
  }

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFiles])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowId]] = {

    val asStoreEntries = sources map toWorkflowStoreEntry
    val returnValue = asStoreEntries map { workflowStore => WorkflowId.fromString(workflowStore.workflowExecutionUuid) }

    // The results from the Future aren't useful, so on completion map it into the precalculated return value instead. Magic!
    sqlDatabase.addWorkflowStoreEntries(asStoreEntries.list.toList) map { _ => returnValue }
  }

  private def fromWorkflowStoreEntry(workflowStoreEntry: WorkflowStoreEntry): WorkflowToStart = {
    val sources = WorkflowSourceFiles(
      workflowStoreEntry.workflowDefinition.toRawString,
      workflowStoreEntry.workflowInputs.toRawString,
      workflowStoreEntry.workflowOptions.toRawString)
    WorkflowToStart(
      WorkflowId.fromString(workflowStoreEntry.workflowExecutionUuid),
      sources,
      fromDbStateStringToStartableState(workflowStoreEntry.workflowState))
  }

  private def toWorkflowStoreEntry(workflowSourceFiles: WorkflowSourceFiles): WorkflowStoreEntry = {
    WorkflowStoreEntry(
      WorkflowId.randomId().toString,
      workflowSourceFiles.wdlSource.toClob,
      workflowSourceFiles.inputsJson.toClob,
      workflowSourceFiles.workflowOptionsJson.toClob,
      WorkflowStoreState.Submitted.toString,
      OffsetDateTime.now.toSystemTimestamp
    )
  }

  private def fromDbStateStringToStartableState(dbStateName: String): StartableState = {
    if (dbStateName.equalsIgnoreCase(WorkflowStoreState.Restartable.toString)) {
      WorkflowStoreState.Restartable
    }
    else {
      WorkflowStoreState.Submitted
    }
  }
}
