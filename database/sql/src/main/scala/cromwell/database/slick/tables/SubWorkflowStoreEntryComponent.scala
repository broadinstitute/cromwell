package cromwell.database.slick.tables

import cromwell.database.sql.tables.SubWorkflowStoreEntry

trait SubWorkflowStoreEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class SubWorkflowStoreEntries(tag: Tag) extends Table[SubWorkflowStoreEntry](tag, "SUB_WORKFLOW_STORE_ENTRY") {
    def subWorkflowStoreEntryId = column[Int]("SUB_WORKFLOW_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def parentWorkflowExecutionUuid = column[String]("PARENT_WORKFLOW_EXECUTION_UUID")

    def callFullyQualifiedName = column[String]("CALL_FULLY_QUALIFIED_NAME")

    def callIndex = column[Int]("CALL_INDEX")

    def callAttempt = column[Int]("CALL_ATTEMPT")

    def subWorkflowExecutionUuid = column[String]("SUB_WORKFLOW_EXECUTION_UUID")

    override def * = (parentWorkflowExecutionUuid, callFullyQualifiedName, callIndex, callAttempt, subWorkflowExecutionUuid, subWorkflowStoreEntryId.?) <> (SubWorkflowStoreEntry.tupled, SubWorkflowStoreEntry.unapply)

    def ucSubWorkflowStoreEntryPweuCfqnJiJa = index("UC_SUB_WORKFLOW_STORE_ENTRY_PWEU_CFQN_CI_CA",
      (parentWorkflowExecutionUuid, callFullyQualifiedName, callIndex, callAttempt), unique = true)

    def ixSubWorkflowStoreEntryPweu = index("IX_SUB_WORKFLOW_STORE_ENTRY_PWEU", parentWorkflowExecutionUuid, unique = false)
  }

  protected val subWorkflowStoreEntries = TableQuery[SubWorkflowStoreEntries]

  val subWorkflowStoreEntryIdsAutoInc = subWorkflowStoreEntries returning subWorkflowStoreEntries.map(_.subWorkflowStoreEntryId)

  val subWorkflowStoreEntriesForParentWorkflowExecutionUuid = Compiled(
    (parentWorkflowExecutionUuid: Rep[String]) => for {
      subWorkflowStoreEntry <- subWorkflowStoreEntries
      if subWorkflowStoreEntry.parentWorkflowExecutionUuid === parentWorkflowExecutionUuid
    } yield subWorkflowStoreEntry
  )

  /**
    * Useful for finding the unique job store for a given job key
    */
  val subWorkflowStoreEntriesForJobKey = Compiled(
    (parentWorkflowExecutionUuid: Rep[String], callFullyQualifiedName: Rep[String], jobIndex: Rep[Int],
     jobAttempt: Rep[Int]) =>
      for {
        subWorkflowStoreEntry <- subWorkflowStoreEntries
        if subWorkflowStoreEntry.parentWorkflowExecutionUuid === parentWorkflowExecutionUuid &&
          subWorkflowStoreEntry.callFullyQualifiedName === callFullyQualifiedName &&
          subWorkflowStoreEntry.callIndex === jobIndex && subWorkflowStoreEntry.callAttempt === jobAttempt
      } yield subWorkflowStoreEntry
  )
}
