package cromwell.database.slick.tables

import java.sql.{Clob, Timestamp}

import cromwell.database.sql.tables.WorkflowStoreEntry

trait WorkflowStoreComponent {

  this: DriverComponent =>

  import driver.api._

  class WorkflowStoreEntries(tag: Tag) extends Table[WorkflowStoreEntry](tag, "WORKFLOW_STORE") {
    def workflowStoreTableId = column[Int]("WORKFLOW_STORE_ID", O.PrimaryKey, O.AutoInc)
    def workflowUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def workflowDefinition = column[Clob]("WORKFLOW_DEFINITION")
    def workflowInputs = column[Clob]("WORKFLOW_INPUTS")
    def workflowOptions = column[Clob]("WORKFLOW_OPTIONS")
    def state = column[String]("WORKFLOW_STATE")
    def submissionTime = column[Timestamp]("SUBMISSION_TIME")

    override def * = (workflowUuid, workflowDefinition, workflowInputs, workflowOptions, state, submissionTime, workflowStoreTableId.?) <>
      (WorkflowStoreEntry.tupled, WorkflowStoreEntry.unapply)

    def uuidIndex = index("WORKFLOW_STORE_UUID_IDX", workflowUuid, unique = true)

    def stateIndex = index("WORKFLOW_STORE_STATE_IDX", state, unique = false)
  }

  protected val workflowStore = TableQuery[WorkflowStoreEntries]

  val workflowStoreAutoInc = workflowStore returning workflowStore.map(_.workflowStoreTableId)

  /**
    * Useful for finding the store entry for a given workflow UUID
    */
  val workflowStoreEntryByWorkflowUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      entry <- workflowStore
      if entry.workflowUuid === workflowExecutionUuid
    } yield entry)

  /**
    * Useful for selecting store entries with a given state.
    */
  val workflowStoreEntriesByState = Compiled(
    (state: Rep[String], limit: ConstColumn[Long]) => {
      val workflowStoreRows = for {
        workflowStoreRow <- workflowStore
        if workflowStoreRow.state === state
      } yield workflowStoreRow
      workflowStoreRows.sortBy(_.submissionTime.asc).take(limit)
    }
  )

  /**
    * Useful for updating state for all entries matching a given UUID
    */
  val workflowStateByUuid = Compiled(
    (uuid: Rep[String]) => for {
      entry <- workflowStore
      if entry.workflowUuid === uuid
    } yield entry.state)

  /**
    * Useful for updating state for all entries matching a given state
    */
  val workflowStoreStateByWorkflowStoreState = Compiled(
    (workflowStoreState: Rep[String]) => for {
      workflowStoreRow <- workflowStore
      if workflowStoreRow.state === workflowStoreState
    } yield workflowStoreRow.state)
}
