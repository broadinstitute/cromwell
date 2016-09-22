package cromwell.database.migration.metadata.table

import cromwell.database.migration.metadata.MetadataCustomSql

class FailureEventTableMigration extends MetadataCustomSql {
  import MetadataCustomSql._

  override def queries: Array[String] = {
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
        |SELECT
        |  WORKFLOW_EXECUTION_UUID,
        |  CONCAT('failures[', fe.FAILURE_EVENT_ID ,']:failure'),
        |  CALL_FQN,
        |  IDX,
        |  ATTEMPT,
        |  fe.EVENT_MESSAGE,
        |  'string',
        |  NOW()
        |FROM FAILURE_EVENT fe
        |  LEFT JOIN EXECUTION e ON fe.EXECUTION_ID = e.EXECUTION_ID
        |  JOIN WORKFLOW_EXECUTION we ON we.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID;
      """.stripMargin,
      s"""
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
        |SELECT
        |  WORKFLOW_EXECUTION_UUID,
        |  CONCAT('failures[', fe.FAILURE_EVENT_ID ,']:timestamp'),
        |  CALL_FQN,
        |  IDX,
        |  ATTEMPT,
        |  DATE_FORMAT(fe.EVENT_TIMESTAMP, '%Y-%m-%dT%T.%f$Offset'),
        |  'string',
        |  NOW()
        |FROM FAILURE_EVENT fe
        |  LEFT JOIN EXECUTION e ON fe.EXECUTION_ID = e.EXECUTION_ID
        |  JOIN WORKFLOW_EXECUTION we ON we.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID;
      """.stripMargin
    )
  }

  override def getConfirmationMessage: String = "Failure Table migration complete."
}
