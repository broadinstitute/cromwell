package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.WorkflowStoreEntry
import slick.profile.RelationalProfile.ColumnOption.Default

trait WorkflowStoreComponent {

  this: DriverComponent =>

  import driver.api._

  class WorkflowStoreEntries(tag: Tag) extends Table[WorkflowStoreEntry](tag, "WORKFLOW_STORE") {
    def workflowStoreTableId = column[Int]("WORKFLOW_STORE_ID", O.PrimaryKey, O.AutoInc)
    def workflowUuid = column[String]("WORKFLOW_UUID")
    def workflowDefinition = column[String]("WORKFLOW_DEFINITION")
    def workflowInputs = column[Option[String]]("WORKFLOW_INPUTS")
    def workflowOptions = column[Option[String]]("WORKFLOW_OPTIONS")
    def state = column[String]("STATE", Default("Submitted"))
    def timestamp = column[Timestamp]("TIMESTAMP")

    override def * = (workflowUuid, workflowDefinition, workflowInputs, workflowOptions, state, timestamp, workflowStoreTableId.?) <>
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
      workflowStoreRows.sortBy(_.timestamp.asc).take(limit)
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
