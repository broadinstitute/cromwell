package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import common.validation.ErrorOr.ErrorOr
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
      // Workflows in Running or Aborting state get their restarted flag set to true on restart  
      sqlDatabase.updateToRestartedIfInState(
        List(WorkflowStoreState.Running.stateName, WorkflowStoreState.Aborting.stateName)
      )
    } else {
      Future.successful(())
    }
  }

  override def aborting(id: WorkflowId)(implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    sqlDatabase.updateWorkflowState(
      id.toString, WorkflowStoreState.Aborting.stateName
    )
  }
  
  override def abortAllRunning()(implicit ec: ExecutionContext): Future[Unit] = {
    sqlDatabase.updateWorkflowsInState(
      List(
        WorkflowStoreState.Running.toString -> WorkflowStoreState.Aborting.toString
      )
    )
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
    import cats.syntax.traverse._
    import cats.instances.list._
    import common.validation.Validation._
    sqlDatabase.queryWorkflowStoreEntries(n, state.stateName, state.restarted, state.afterFetchedState.stateName) map {
      // .get on purpose here to fail the future if something went wrong
      _.toList.traverse(fromWorkflowStoreEntry).toTry.get
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

  private def fromWorkflowStoreEntry(workflowStoreEntry: WorkflowStoreEntry): ErrorOr[WorkflowToStart] = {
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

    fromDbStateStringToStartableState(workflowStoreEntry.workflowState, workflowStoreEntry.restarted) map { startableState =>
      WorkflowToStart(
        WorkflowId.fromString(workflowStoreEntry.workflowExecutionUuid),
        sources,
        startableState)
    }
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
      restarted = false,
      submissionTime = OffsetDateTime.now.toSystemTimestamp,
      importsZip = workflowSourceFiles.importsZipFileOption.toBlobOption
    )
  }

  private def fromDbStateStringToStartableState(dbStateName: String, restarted: Boolean): ErrorOr[StartableState] = {
    import cats.syntax.validated._

    WorkflowStoreState.fromNameAndRestarted(dbStateName, restarted) match {
      case Some(startable: StartableState) => startable.validNel
        // Both of those should be impossible unless someone is messing with the database
      case Some(nonStartable) => s"$nonStartable is not a startable state".invalidNel
      case None => s"Cannot find a known state for $dbStateName".invalidNel
    }
  }
}
