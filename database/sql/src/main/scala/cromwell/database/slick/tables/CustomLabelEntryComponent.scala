package cromwell.database.slick.tables

import cromwell.database.sql.tables.CustomLabelEntry
import slick.model.ForeignKeyAction.Cascade

trait CustomLabelEntryComponent {

  this: DriverComponent with WorkflowMetadataSummaryEntryComponent =>

  import driver.api._

  class CustomLabelEntries(tag: Tag)
    extends Table[CustomLabelEntry](tag, "CUSTOM_LABEL_ENTRY") {
    def customLabelEntryId = column[Long]("CUSTOM_LABEL_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def customLabelKey = column[String]("CUSTOM_LABEL_KEY", O.Length(63))

    def customLabelValue = column[String]("CUSTOM_LABEL_VALUE", O.Length(63))

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(100))

    override def * = (customLabelKey, customLabelValue, workflowExecutionUuid,
    customLabelEntryId.?) <> (CustomLabelEntry.tupled, CustomLabelEntry.unapply)

    def fkCustomLabelEntryWorkflowExecutionUuid = foreignKey("FK_CUSTOM_LABEL_ENTRY_WORKFLOW_EXECUTION_UUID",
      workflowExecutionUuid, workflowMetadataSummaryEntries)(_.workflowExecutionUuid, onDelete = Cascade)

    def ucCustomLabelEntryClkClvWeu = index("UC_CUSTOM_LABEL_ENTRY_CLK_CLV_WEU",
      (customLabelKey, customLabelValue, workflowExecutionUuid), unique = true)
}

  val customLabelEntries = TableQuery[CustomLabelEntries]

  val customLabelEntryIdsAutoInc = customLabelEntries returning
    customLabelEntries.map(_.customLabelEntryId)

  val existsWorkflowIdLabelKeyAndValue = Compiled(
    (workflowUuid: Rep[String], labelKey: Rep[String], labelValue: Rep[String]) => (for {
      customLabelEntry <- customLabelEntries
      if customLabelEntry.workflowExecutionUuid === workflowUuid &&
        customLabelEntry.customLabelKey == labelKey &&
        customLabelEntry.customLabelValue == labelValue
    } yield ()).exists
  )

  def existsWorkflowIdLabelKeyAndValue(workflowId: Rep[String],
                                       labelKey: Rep[String],
                                       labelValue: Rep[String]): Rep[Boolean] = {
    customLabelEntries.filter(customLabelEntry =>
      customLabelEntry.workflowExecutionUuid === workflowId &&
        customLabelEntry.customLabelKey === labelKey &&
        customLabelEntry.customLabelValue === labelValue
    ).exists
  }
}
