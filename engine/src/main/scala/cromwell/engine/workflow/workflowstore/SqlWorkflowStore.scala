package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.WorkflowStoreSqlDatabase
import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState
import eu.timepit.refined.api.Refined
import eu.timepit.refined.collection._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

case class SqlWorkflowStore(sqlDatabase: WorkflowStoreSqlDatabase) extends WorkflowStore {
  override def initialize(implicit ec: ExecutionContext): Future[Unit] = {
    if (ConfigFactory.load().as[Option[Boolean]]("system.workflow-restart").getOrElse(true)) {
      sqlDatabase.updateWorkflowsInState(
        List(
          WorkflowStoreState.Running.toString -> WorkflowStoreState.RestartableRunning.toString,
          WorkflowStoreState.Aborting.toString -> WorkflowStoreState.RestartableAborting.toString
        )
      )
    } else {
      Future.successful(())
    }
  }

  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    sqlDatabase.updateWorkflowState(
      id.toString, WorkflowStoreState.Aborting.toString
    ).map(_ > 0)
  }

  override def stats(implicit ec: ExecutionContext): Future[Map[String, Int]] = sqlDatabase.stats

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = {
    sqlDatabase.removeWorkflowStoreEntry(id.toString).map(_ > 0) // i.e. did anything get deleted
  }

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  override def fetchRunnableWorkflows(n: Int, state: StartableState)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    sqlDatabase.queryWorkflowStoreEntries(n, state.toString, state.afterFetchedState.toString) map {
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
      workflowSource = workflowStoreEntry.workflowDefinition.toRawString,
      workflowType = workflowStoreEntry.workflowType,
      workflowTypeVersion = workflowStoreEntry.workflowTypeVersion,
      inputsJson = workflowStoreEntry.workflowInputs.toRawString,
      workflowOptionsJson = workflowStoreEntry.workflowOptions.toRawString,
      labelsJson = workflowStoreEntry.customLabels.toRawString,
      importsFile = workflowStoreEntry.importsZip.toBytesOption,
      warnings = Vector.empty
    )
    WorkflowToStart(
      WorkflowId.fromString(workflowStoreEntry.workflowExecutionUuid),
      sources,
      fromDbStateStringToStartableState(workflowStoreEntry.workflowState))
  }

  private def toWorkflowStoreEntry(workflowSourceFiles: WorkflowSourceFilesCollection): WorkflowStoreEntry = {
    import eu.timepit.refined._
    val nonEmptyJsonString: String Refined NonEmpty  = refineMV[NonEmpty]("{}")

    WorkflowStoreEntry(
      workflowExecutionUuid = WorkflowId.randomId().toString,
      workflowDefinition = workflowSourceFiles.workflowSource.toClobOption,
      workflowType = workflowSourceFiles.workflowType,
      workflowTypeVersion = workflowSourceFiles.workflowTypeVersion,
      workflowInputs = workflowSourceFiles.inputsJson.toClobOption,
      workflowOptions = workflowSourceFiles.workflowOptionsJson.toClobOption,
      customLabels = workflowSourceFiles.labelsJson.toClob(default = nonEmptyJsonString),
      workflowState = WorkflowStoreState.Submitted.toString,
      submissionTime = OffsetDateTime.now.toSystemTimestamp,
      importsZip = workflowSourceFiles.importsZipFileOption.toBlobOption
    )
  }

  private def fromDbStateStringToStartableState(dbStateName: String): StartableState = {
    if (dbStateName.equalsIgnoreCase(WorkflowStoreState.RestartableRunning.toString)) {
      WorkflowStoreState.RestartableRunning
    } else if(dbStateName.equalsIgnoreCase(WorkflowStoreState.RestartableAborting.toString)){
      WorkflowStoreState.RestartableAborting
    } else {
      WorkflowStoreState.Submitted
    }
  }
}
