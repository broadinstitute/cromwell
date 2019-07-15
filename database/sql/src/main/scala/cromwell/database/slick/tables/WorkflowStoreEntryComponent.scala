package cromwell.database.slick.tables

import java.sql.Timestamp
import javax.sql.rowset.serial.{SerialBlob, SerialClob}
import cromwell.database.sql.tables.WorkflowStoreEntry

trait WorkflowStoreEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class WorkflowStoreEntries(tag: Tag) extends Table[WorkflowStoreEntry](tag, "WORKFLOW_STORE_ENTRY") {
    def workflowStoreEntryId = column[Int]("WORKFLOW_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(255))

    def workflowRoot = column[Option[String]]("WORKFLOW_ROOT", O.Length(100))

    def workflowType = column[Option[String]]("WORKFLOW_TYPE", O.Length(30))

    def workflowTypeVersion = column[Option[String]]("WORKFLOW_TYPE_VERSION", O.Length(255))

    def workflowDefinition = column[Option[SerialClob]]("WORKFLOW_DEFINITION")

    def workflowUrl = column[Option[String]]("WORKFLOW_URL", O.Length(2000))

    def workflowInputs = column[Option[SerialClob]]("WORKFLOW_INPUTS")

    def workflowOptions = column[Option[SerialClob]]("WORKFLOW_OPTIONS")

    def customLabels = column[SerialClob]("CUSTOM_LABELS")

    def workflowState = column[String]("WORKFLOW_STATE", O.Length(20))

    def submissionTime = column[Timestamp]("SUBMISSION_TIME")

    def importsZip = column[Option[SerialBlob]]("IMPORTS_ZIP")

    def cromwellId = column[Option[String]]("CROMWELL_ID", O.Length(100))

    def heartbeatTimestamp = column[Option[Timestamp]]("HEARTBEAT_TIMESTAMP")

    def hogGroup = column[Option[String]]("HOG_GROUP", O.Length(100))

    override def * = (workflowExecutionUuid, workflowDefinition, workflowUrl, workflowRoot, workflowType, workflowTypeVersion, workflowInputs, workflowOptions, workflowState,
      submissionTime, importsZip, customLabels, cromwellId, heartbeatTimestamp, hogGroup, workflowStoreEntryId.?) <> ((WorkflowStoreEntry.apply _).tupled, WorkflowStoreEntry.unapply)

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
    (limit: ConstColumn[Long],
     heartbeatTimestampTimedOut: ConstColumn[Timestamp],
     excludeWorkflowState: Rep[String]
    ) => {
      val query = for {
        row <- workflowStoreEntries
        /*
        This looks for:

        1) Workflows with no heartbeat (newly submitted or from a cleanly shut down Cromwell).
        2) Workflows with old heartbeats, presumably abandoned by a defunct Cromwell.

        Workflows are taken by submission time, oldest first. This is a "query for update", meaning rows are
        locked such that readers are blocked since we will do an update subsequent to this select in the same
        transaction that we know will impact those readers.
         */
        if (row.heartbeatTimestamp.isEmpty || row.heartbeatTimestamp < heartbeatTimestampTimedOut) &&
          (row.workflowState =!= excludeWorkflowState)
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
    (workflowState: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowState === workflowState
    } yield workflowStoreEntry.workflowState
  )

  /**
    * Useful for updating a given workflow to a new state
    */
  val workflowStateForWorkflowExecutionUUid = Compiled(
    (workflowId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowExecutionUuid === workflowId
    } yield workflowStoreEntry.workflowState
  )

  /**
    * Useful for updating a given workflow to a 'Submitted' state when it's currently 'On Hold'
    */
  val workflowStateForWorkflowExecutionUUidAndWorkflowState = Compiled(
    (workflowId: Rep[String], workflowState: Rep[String]) => {
      for {
        workflowStoreEntry <- workflowStoreEntries
        if workflowStoreEntry.workflowExecutionUuid === workflowId
        if workflowStoreEntry.workflowState === workflowState
      } yield workflowStoreEntry.workflowState
    }
  )

  /**
    * Useful for deleting a given workflow to a 'Submitted' state when it's currently 'On Hold' or 'Submitted'
    */
  val workflowStoreEntryForWorkflowExecutionUUidAndWorkflowStates = Compiled(
    (workflowId: Rep[String],
     workflowStateOr1: Rep[String],
     workflowStateOr2: Rep[String]
    ) => {
      for {
        workflowStoreEntry <- workflowStoreEntries
        if workflowStoreEntry.workflowExecutionUuid === workflowId
        if workflowStoreEntry.workflowState === workflowStateOr1 ||
          workflowStoreEntry.workflowState === workflowStateOr2
      } yield workflowStoreEntry
    }
  )

  val findWorkflowsWithAbortRequested = Compiled(
    (cromwellId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.workflowState === "Aborting" && workflowStoreEntry.cromwellId === cromwellId
    } yield workflowStoreEntry.workflowExecutionUuid
  )

  val findWorkflows = Compiled(
    (cromwellId: Rep[String]) => for {
      workflowStoreEntry <- workflowStoreEntries
      if workflowStoreEntry.cromwellId === cromwellId
    } yield workflowStoreEntry.workflowExecutionUuid
  )
}
