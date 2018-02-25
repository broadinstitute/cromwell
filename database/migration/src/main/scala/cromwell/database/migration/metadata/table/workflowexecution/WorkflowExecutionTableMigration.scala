package cromwell.database.migration.metadata.table.workflowexecution

import cromwell.database.migration.metadata.MetadataCustomSql
import MetadataCustomSql._

class WorkflowExecutionTableMigration extends MetadataCustomSql {

  override def queries: Array[String] = {
    Array(
    s"""
      |INSERT INTO METADATA_JOURNAL (
      |WORKFLOW_EXECUTION_UUID,
      |METADATA_KEY,
      |METADATA_VALUE,
      |METADATA_VALUE_TYPE,
      |METADATA_TIMESTAMP
      |)
      |SELECT
      |  WORKFLOW_EXECUTION_UUID,
      |  'submission',
      |  DATE_FORMAT(START_DT, '%Y-%m-%dT%T.%f$Offset'),
      |  'string',
      |  NOW()
      |FROM WORKFLOW_EXECUTION
      |WHERE START_DT IS NOT NULL;""".stripMargin
    ,
    s"""
      |INSERT INTO METADATA_JOURNAL (
      |WORKFLOW_EXECUTION_UUID,
      |METADATA_KEY,
      |METADATA_VALUE,
      |METADATA_VALUE_TYPE,
      |METADATA_TIMESTAMP
      |)
      |SELECT
      |  WORKFLOW_EXECUTION_UUID,
      |  'start',
      |  DATE_FORMAT(START_DT, '%Y-%m-%dT%T.%f$Offset'),
      |  'string',
      |  NOW()
      |FROM WORKFLOW_EXECUTION
      |WHERE START_DT IS NOT NULL;""".stripMargin
    ,
    s"""
      |INSERT INTO METADATA_JOURNAL (
      |WORKFLOW_EXECUTION_UUID,
      |METADATA_KEY,
      |METADATA_VALUE,
      |METADATA_VALUE_TYPE,
      |METADATA_TIMESTAMP
      |)
      |SELECT
      |  WORKFLOW_EXECUTION_UUID,
      |  'end',
      |  DATE_FORMAT(END_DT, '%Y-%m-%dT%T.%f$Offset'),
      |  'string',
      |  NOW()
      |FROM WORKFLOW_EXECUTION
      |WHERE END_DT IS NOT NULL;""".stripMargin
    ,
    s"""
      |INSERT INTO METADATA_JOURNAL (
      |WORKFLOW_EXECUTION_UUID,
      |METADATA_KEY,
      |METADATA_VALUE,
      |METADATA_VALUE_TYPE,
      |METADATA_TIMESTAMP
      |)
      |SELECT
      |  WORKFLOW_EXECUTION_UUID,
      |  'status',
      |  STATUS,
      |  'string',
      |  NOW()
      |FROM WORKFLOW_EXECUTION;""".stripMargin
    ,
    """
      |INSERT INTO METADATA_JOURNAL (
      |WORKFLOW_EXECUTION_UUID,
      |METADATA_KEY,
      |METADATA_VALUE,
      |METADATA_VALUE_TYPE,
      |METADATA_TIMESTAMP
      |)
      |SELECT
      |  WORKFLOW_EXECUTION_UUID,
      |  'workflowName',
      |  WORKFLOW_NAME,
      |  'string',
      |  NOW()
      |FROM WORKFLOW_EXECUTION;""".stripMargin
    ,
    """
      |INSERT INTO METADATA_JOURNAL (
      |WORKFLOW_EXECUTION_UUID,
      |METADATA_KEY,
      |METADATA_VALUE,
      |METADATA_VALUE_TYPE,
      |METADATA_TIMESTAMP
      |)
      |SELECT
      |  WORKFLOW_EXECUTION_UUID,
      |  'outputs',
      |  NULL,
      |  NULL,
      |  '1900-01-01 0.000000'
      |FROM WORKFLOW_EXECUTION;""".stripMargin
    )
  }

  override def getConfirmationMessage: String = "Workflow Execution Table migration complete."
}
