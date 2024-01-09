package cromwell.database.slick.tables

import cromwell.database.sql.tables.SubWorkflowStoreEntry
import slick.model.ForeignKeyAction.Cascade

trait SubWorkflowStoreEntryComponent {

  this: DriverComponent with WorkflowStoreEntryComponent =>

  import driver.api._

  class SubWorkflowStoreEntries(tag: Tag) extends Table[SubWorkflowStoreEntry](tag, "SUB_WORKFLOW_STORE_ENTRY") {
    def subWorkflowStoreEntryId = column[Long]("SUB_WORKFLOW_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def rootWorkflowId = column[Long]("ROOT_WORKFLOW_ID")

    def parentWorkflowExecutionUuid = column[String]("PARENT_WORKFLOW_EXECUTION_UUID", O.Length(255))

    def callFullyQualifiedName = column[String]("CALL_FULLY_QUALIFIED_NAME", O.Length(255))

    def callIndex = column[Int]("CALL_INDEX")

    def callAttempt = column[Int]("CALL_ATTEMPT")

    def subWorkflowExecutionUuid = column[String]("SUB_WORKFLOW_EXECUTION_UUID", O.Length(255))

    override def * = (rootWorkflowId.?,
                      parentWorkflowExecutionUuid,
                      callFullyQualifiedName,
                      callIndex,
                      callAttempt,
                      subWorkflowExecutionUuid,
                      subWorkflowStoreEntryId.?
    ) <> (SubWorkflowStoreEntry.tupled, SubWorkflowStoreEntry.unapply)

    def ucSubWorkflowStoreEntryPweuCfqnCiCa = index(
      "UC_SUB_WORKFLOW_STORE_ENTRY_PWEU_CFQN_CI_CA",
      (parentWorkflowExecutionUuid, callFullyQualifiedName, callIndex, callAttempt),
      unique = true
    )

    def fkSubWorkflowStoreEntryRootWorkflowId = foreignKey("FK_SUB_WORKFLOW_STORE_ENTRY_ROOT_WORKFLOW_ID",
                                                           rootWorkflowId,
                                                           workflowStoreEntries
    )(_.workflowStoreEntryId, onDelete = Cascade)

    def ixSubWorkflowStoreEntryPweu =
      index("IX_SUB_WORKFLOW_STORE_ENTRY_PWEU", parentWorkflowExecutionUuid, unique = false)
  }

  protected val subWorkflowStoreEntries = TableQuery[SubWorkflowStoreEntries]

  val subWorkflowStoreEntryIdsAutoInc =
    subWorkflowStoreEntries returning subWorkflowStoreEntries.map(_.subWorkflowStoreEntryId)

  val subWorkflowStoreEntriesForRootWorkflowId = Compiled((rootWorkflowId: Rep[Long]) =>
    for {
      subWorkflowStoreEntry <- subWorkflowStoreEntries
      if subWorkflowStoreEntry.rootWorkflowId === rootWorkflowId
    } yield subWorkflowStoreEntry
  )

  /**
    * Useful for finding the unique sub workflow entry for a given job key
    */
  val subWorkflowStoreEntriesForJobKey = Compiled(
    (parentWorkflowExecutionUuid: Rep[String],
     callFullyQualifiedName: Rep[String],
     jobIndex: Rep[Int],
     jobAttempt: Rep[Int]
    ) =>
      for {
        subWorkflowStoreEntry <- subWorkflowStoreEntries
        if subWorkflowStoreEntry.parentWorkflowExecutionUuid === parentWorkflowExecutionUuid &&
          subWorkflowStoreEntry.callFullyQualifiedName === callFullyQualifiedName &&
          subWorkflowStoreEntry.callIndex === jobIndex && subWorkflowStoreEntry.callAttempt === jobAttempt
      } yield subWorkflowStoreEntry
  )
}
