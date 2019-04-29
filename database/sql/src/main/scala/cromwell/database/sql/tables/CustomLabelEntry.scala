package cromwell.database.sql.tables

case class CustomLabelEntry
(
  customLabelKey: String,
  customLabelValue: String,
  workflowExecutionUuid: String,
  customLabelEntryId: Option[Long] = None
) {
  def toUpdateColumns = (customLabelKey, customLabelValue, workflowExecutionUuid)
}
