package cromwell.database.migration.metadata.table.executionevent

import cromwell.database.migration.metadata.MetadataCustomSql

class ExecutionEventTableDescriptionMigration extends MetadataCustomSql {

  override def queries: Array[String] =
    Array(
      """
        |INSERT INTO METADATA_JOURNAL (
        |  WORKFLOW_EXECUTION_UUID,
        |  METADATA_KEY,
        |  CALL_FQN,
        |  JOB_SCATTER_INDEX,
        |  JOB_RETRY_ATTEMPT,
        |  METADATA_VALUE,
        |  METADATA_VALUE_TYPE,
        |  METADATA_TIMESTAMP
        |)
        |  SELECT
        |    WORKFLOW_EXECUTION_UUID,
        |    CONCAT('executionEvents[', ev.EVENT_ID ,']:description'),
        |    CALL_FQN,
        |    IDX,
        |    ATTEMPT,
        |    ev.DESCRIPTION,
        |    'string',
        |    NOW()
        |  FROM EXECUTION_EVENT ev
        |    LEFT JOIN EXECUTION e ON ev.EXECUTION_ID = e.EXECUTION_ID
        |    JOIN WORKFLOW_EXECUTION we ON we.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID;
      """.stripMargin
    )

  override def getConfirmationMessage: String = "Execution Event Table (Description field) migration complete."
}
