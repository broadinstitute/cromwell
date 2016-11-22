package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime
import javax.sql.rowset.serial.SerialBlob

import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

import scala.concurrent.{ExecutionContext, Future}

case class SqlWorkflowStore(sqlDatabase: WorkflowStoreSqlDatabase) extends WorkflowStore {
  override def initialize(implicit ec: ExecutionContext): Future[Unit] = {
          if (ConfigFactory.load().as[Option[Boolean]]("system.workflow-restart").getOrElse(true)) {
            sqlDatabase.updateWorkflowState(
                WorkflowStoreState.Running.toString,
                WorkflowStoreState.Restartable.toString)
          } else {
            Future.successful(())
          }
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
  override def add(sources: NonEmptyList[WorkflowSourceFilesCollection])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowId]] = {

    val asStoreEntries = sources map toWorkflowStoreEntry
    val returnValue = asStoreEntries map { workflowStore => WorkflowId.fromString(workflowStore.workflowExecutionUuid) }

    // The results from the Future aren't useful, so on completion map it into the precalculated return value instead. Magic!
    sqlDatabase.addWorkflowStoreEntries(asStoreEntries.toList) map { _ => returnValue }
  }

  private def fromWorkflowStoreEntry(workflowStoreEntry: WorkflowStoreEntry): WorkflowToStart = {
    val sources = WorkflowSourceFilesCollection(
      workflowStoreEntry.workflowDefinition.toRawString,
      workflowStoreEntry.workflowInputs.toRawString,
      workflowStoreEntry.workflowOptions.toRawString,
      workflowStoreEntry.importsZipFile.map(b => b.getBytes(1, b.length.asInstanceOf[Int]))
    )
    WorkflowToStart(
      WorkflowId.fromString(workflowStoreEntry.workflowExecutionUuid),
      sources,
      fromDbStateStringToStartableState(workflowStoreEntry.workflowState))
  }

  private def toWorkflowStoreEntry(workflowSourceFiles: WorkflowSourceFilesCollection): WorkflowStoreEntry = {
    WorkflowStoreEntry(
      WorkflowId.randomId().toString,
      workflowSourceFiles.wdlSource.toClob,
      workflowSourceFiles.inputsJson.toClob,
      workflowSourceFiles.workflowOptionsJson.toClob,
      WorkflowStoreState.Submitted.toString,
      OffsetDateTime.now.toSystemTimestamp,
      workflowSourceFiles.importsZipFileOption.map(new SerialBlob(_))
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
