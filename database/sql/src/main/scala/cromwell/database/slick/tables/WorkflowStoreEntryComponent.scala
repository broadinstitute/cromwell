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

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255))

    def workflowRoot = column[Option[String]]("WORKFLOW_ROOT", O.Length(100))

    def workflowType = column[Option[String]]("WORKFLOW_TYPE", O.Length(30))

    def workflowTypeVersion = column[Option[String]]("WORKFLOW_TYPE_VERSION", O.Length(255))

    def workflowDefinition = column[Option[Clob]]("WORKFLOW_DEFINITION")

    def workflowUrl = column[Option[String]]("WORKFLOW_URL", O.Length(2000))

    def workflowInputs = column[Option[Clob]]("WORKFLOW_INPUTS")

    def workflowOptions = column[Option[Clob]]("WORKFLOW_OPTIONS")

    def customLabels = column[Clob]("CUSTOM_LABELS")

    def workflowState = column[WorkflowStoreState]("WORKFLOW_STATE", O.Length(20))

    def submissionTime = column[Timestamp]("SUBMISSION_TIME")

    def importsZip = column[Option[Blob]]("IMPORTS_ZIP")

    def cromwellId = column[Option[String]]("CROMWELL_ID", O.Length(100))

    def heartbeatTimestamp = column[Option[Timestamp]]("HEARTBEAT_TIMESTAMP")

    override def * = (workflowExecutionUuid, workflowDefinition, workflowUrl, workflowRoot, workflowType, workflowTypeVersion, workflowInputs, workflowOptions, workflowState,
      submissionTime, importsZip, customLabels, cromwellId, heartbeatTimestamp, workflowStoreEntryId.?) <> ((WorkflowStoreEntry.apply _).tupled, WorkflowStoreEntry.unapply)

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

  val heartbeatForWorkflowStoreEntry = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowStoreEntry.heartbeatTimestamp
  )

  /**
    * Returns up to "limit" startable workflows, sorted by submission time.
    */
  val fetchStartableWorkflows = Compiled(
    (limit: ConstColumn[Long], heartbeatThreshold: ConstColumn[Timestamp]) => {
      val query = for {
        row <- workflowStoreEntries
        // This looks for:
        //
        // 1) Workflows with no heartbeat (newly submitted or from a cleanly shut down Cromwell).
        // 2) Workflows with old heartbeats, presumably abandoned by a defunct Cromwell.
        //
        // Workflows are taken by submission time, oldest first. This is a "query for update", meaning rows are
        // locked such that readers are blocked since we will do an update subsequent to this select in the same
        // transaction that we know will impact those readers.
        //
        // The current code only writes heartbeats on initial pickup; this needs to be fixed by instrumenting workflow
        // heartbeats into Cromwell appropriately.  But in the meantime any other Cromwell's workflows will appear
        // to be abandoned if they were picked up before `heartbeatThreshold`.
        if (row.heartbeatTimestamp.isEmpty || row.heartbeatTimestamp < heartbeatThreshold) && !(row.workflowState === WorkflowStoreState.OnHold)
      } yield row
      query.forUpdate.sortBy(_.submissionTime.asc).take(limit)
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
    * Useful for updating the relevant fields of a workflow store entry when a workflow is picked up for processing.
    */
  val workflowStoreFieldsForPickup = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      row <- workflowStoreEntries
      if row.workflowExecutionUuid === workflowExecutionUuid
    } yield (row.workflowState, row.cromwellId, row.heartbeatTimestamp)
  )

  /**
    * Useful for clearing out cromwellId and heartbeatTimestamp on an orderly Cromwell shutdown.
    */
  val releaseWorkflowStoreEntries = Compiled(
    (cromwellId: Rep[String]) => for {
      row <- workflowStoreEntries
      if row.cromwellId === cromwellId
    } yield (row.cromwellId, row.heartbeatTimestamp)
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
    * Useful for clearing the heartbeat timestamp on server restart.
    */
  val heartbeatTimestamp = Compiled(
    (cromwellId: Rep[Option[String]]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if (workflowStoreEntry.cromwellId === cromwellId) || (workflowStoreEntry.cromwellId.isEmpty && cromwellId.isEmpty)
    } yield workflowStoreEntry.heartbeatTimestamp
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
    * Useful for updating a given workflow to a 'Submitted' state when it's currently 'On Hold'
    */
  val workflowForIdAndOnHold = Compiled(
    (workflowId: Rep[String]) => {
       for {
        workflowStoreEntry <- workflowStoreEntries
        if workflowStoreEntry.workflowExecutionUuid === workflowId
        if workflowStoreEntry.workflowState === WorkflowStoreState.OnHold
      } yield workflowStoreEntry.workflowState
    }
  )

  /**
    * Useful for checking if the heartbeat timestamp was cleared for a given workflow id.
    */
  val heartbeatClearedForWorkflowId = Compiled(
    (workflowId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowId
    } yield workflowStoreEntry.heartbeatTimestamp.isEmpty
  )
}
