package cromwell.database.slick.tables

import java.sql.{Blob, Clob, Timestamp}

import cromwell.database.sql.tables.WorkflowStoreEntry
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState

trait WorkflowStoreEntryComponent {

  this: DriverComponent =>

  import driver.api._

  object WorkflowStoreEntries {
    implicit val workflowStoreStateMapper = MappedColumnType.base[WorkflowStoreState, String](
      e => e.toString,
      s => WorkflowStoreState.withName(s)
    )
  }

  import WorkflowStoreEntries._
  
  class WorkflowStoreEntries(tag: Tag) extends Table[WorkflowStoreEntry](tag, "WORKFLOW_STORE_ENTRY") {
    def workflowStoreEntryId = column[Int]("WORKFLOW_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def workflowType = column[Option[String]]("WORKFLOW_TYPE", O.Length(30))

    def workflowTypeVersion = column[Option[String]]("WORKFLOW_TYPE_VERSION")

    def workflowDefinition = column[Option[Clob]]("WORKFLOW_DEFINITION")

    def workflowInputs = column[Option[Clob]]("WORKFLOW_INPUTS")

    def workflowOptions = column[Option[Clob]]("WORKFLOW_OPTIONS")

    def customLabels = column[Clob]("CUSTOM_LABELS")

    def workflowState = column[WorkflowStoreState]("WORKFLOW_STATE", O.Length(20))

    def restarted = column[Boolean]("RESTARTED")

    def submissionTime = column[Timestamp]("SUBMISSION_TIME")

    def importsZip = column[Option[Blob]]("IMPORTS_ZIP")

    override def * = (workflowExecutionUuid, workflowDefinition, workflowType, workflowTypeVersion, workflowInputs, workflowOptions, workflowState,
      restarted, submissionTime, importsZip, customLabels, workflowStoreEntryId.?) <> ((WorkflowStoreEntry.apply _).tupled, WorkflowStoreEntry.unapply)

    def ucWorkflowStoreEntryWeu = index("UC_WORKFLOW_STORE_ENTRY_WEU", workflowExecutionUuid, unique = true)

    def ixWorkflowStoreEntryWs = index("IX_WORKFLOW_STORE_ENTRY_WS", workflowState, unique = false)
  }

  protected val workflowStoreEntries = TableQuery[WorkflowStoreEntries]

  val workflowStoreEntryIdsAutoInc = workflowStoreEntries returning workflowStoreEntries.map(_.workflowStoreEntryId)

  /**
    * Useful for finding the workflow store for a given workflow execution UUID
    */
  val workflowStoreEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowStoreEntry
  )

  /**
    * Returns up to "limit" startable workflows, sorted by submission time.
    */
  val fetchStartableWorkflows = Compiled(
    (limit: ConstColumn[Long]) => {
      val query = for {
        workflowStoreEntryRow <- workflowStoreEntries
        if (workflowStoreEntryRow.workflowState === WorkflowStoreState.Aborting && workflowStoreEntryRow.restarted === true) ||
           (workflowStoreEntryRow.workflowState === WorkflowStoreState.Running && workflowStoreEntryRow.restarted === true) ||
           (workflowStoreEntryRow.workflowState === WorkflowStoreState.Submitted && workflowStoreEntryRow.restarted === false)
      } yield workflowStoreEntryRow
      query.sortBy(_.submissionTime.asc).take(limit)
    }
  )

  /**
    * Useful for counting workflows in a given state.
    */
  val workflowStoreStats = Compiled(
    for {
      (state, entry) <- workflowStoreEntries groupBy (_.workflowState)
    } yield state -> entry.size
  )

  /**
    * Useful for updating state for all entries matching a given UUID
    */
  val workflowStateAndRestartedForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowStoreEntry.workflowState, workflowStoreEntry.restarted)
  )

  /**
    * Useful for updating state for all entries matching a given state
    */
  val workflowStateForWorkflowState = Compiled(
    (workflowState: Rep[WorkflowStoreState]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowState === workflowState
    } yield workflowStoreEntry.workflowState
  )

  /**
    * Useful for updating restarted flags on server restart.
    */
  val restartedFlagForRunningAndAborting = Compiled(
    for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowState === WorkflowStoreState.Running || workflowStoreEntry.workflowState === WorkflowStoreState.Aborting
    } yield workflowStoreEntry.restarted
  )

  /**
    * Useful for updating a given workflow to a new state
    */
  val workflowStateForId = Compiled(
    (workflowId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowId
    } yield workflowStoreEntry.workflowState
  )

  /**
    * Useful for updating the restarted flag for a given workflow
    */
  val workflowRestartedForId = Compiled(
    (workflowId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowId
    } yield workflowStoreEntry.restarted
  )
}
