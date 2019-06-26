package cromwell.database.slick.tables

import cromwell.database.sql.tables.CustomLabelEntry
import shapeless.syntax.std.tuple._
import slick.model.ForeignKeyAction.Cascade

trait CustomLabelEntryComponent {

  this: DriverComponent with WorkflowMetadataSummaryEntryComponent =>

  import driver.api.TupleMethods._
  import driver.api._

  class CustomLabelEntries(tag: Tag)
    extends Table[CustomLabelEntry](tag, "CUSTOM_LABEL_ENTRY") {
    def customLabelEntryId = column[Long]("CUSTOM_LABEL_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def customLabelKey = column[String]("CUSTOM_LABEL_KEY", O.Length(255))

    def customLabelValue = column[String]("CUSTOM_LABEL_VALUE", O.Length(255))

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(100))

    def baseProjection = (customLabelKey, customLabelValue, workflowExecutionUuid)

    override def * = baseProjection ~ customLabelEntryId.? <> (CustomLabelEntry.tupled, CustomLabelEntry.unapply)

    def forUpdate = baseProjection.shaped <> (
      tuple => CustomLabelEntry.tupled(tuple :+ None),
      CustomLabelEntry.unapply(_: CustomLabelEntry).map(_.reverse.tail.reverse)
    )

    def fkCustomLabelEntryWorkflowExecutionUuid = foreignKey("FK_CUSTOM_LABEL_ENTRY_WORKFLOW_EXECUTION_UUID",
      workflowExecutionUuid, workflowMetadataSummaryEntries)(_.workflowExecutionUuid, onDelete = Cascade)

    def ucCustomLabelEntryClkWeu = index("UC_CUSTOM_LABEL_ENTRY_CLK_WEU",
      (customLabelKey, workflowExecutionUuid), unique = true)

    def ixCustomLabelEntryClkClv = index("IX_CUSTOM_LABEL_ENTRY_CLK_CLV", (customLabelKey, customLabelValue), unique = false)
}

  val customLabelEntries = TableQuery[CustomLabelEntries]

  val customLabelEntryIdsAutoInc = customLabelEntries returning
    customLabelEntries.map(_.customLabelEntryId)

  val customLabelEntriesForWorkflowExecutionUuidAndLabelKey = Compiled(
    (workflowExecutionUuid: Rep[String], labelKey: Rep[String]) => for {
      customLabelEntry <- customLabelEntries
      if customLabelEntry.workflowExecutionUuid === workflowExecutionUuid &&
        customLabelEntry.customLabelKey === labelKey
    } yield customLabelEntry.forUpdate)

  def existsWorkflowIdLabelKeyAndValue(workflowId: Rep[String],
                                       labelKey: Rep[String],
                                       labelValue: Rep[String]): Rep[Boolean] = {
    customLabelEntries.filter(customLabelEntry =>
      customLabelEntry.workflowExecutionUuid === workflowId &&
        customLabelEntry.customLabelKey === labelKey &&
        customLabelEntry.customLabelValue === labelValue
    ).exists
  }

  val labelsForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      customLabelEntry <- customLabelEntries
      if customLabelEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield (customLabelEntry.customLabelKey, customLabelEntry.customLabelValue)
  )
}
