package cromwell.database.migration.metadata.table.symbol

import java.sql.{PreparedStatement, Timestamp, Types}
import java.time.format.DateTimeFormatter
import java.time.{OffsetDateTime, ZoneId, ZoneOffset}

import org.slf4j.LoggerFactory
import wdl4s.values.{WdlBoolean, WdlFloat, WdlInteger, WdlValue}

object MetadataStatement {
  val WorkflowIdIdx = 1
  val KeyIdx = 2
  val CallFqnIdx = 3
  val CallIndexIdx = 4
  val CallAttemptIdx = 5
  val ValueIdx = 6
  val TimestampIdx = 7
  val ValueTypeIdx = 8

  val InsertSql =
    """
      |INSERT INTO METADATA_JOURNAL
      |(WORKFLOW_EXECUTION_UUID, METADATA_KEY, CALL_FQN, JOB_SCATTER_INDEX, JOB_RETRY_ATTEMPT, METADATA_VALUE, METADATA_TIMESTAMP, METADATA_VALUE_TYPE)
      |VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """.stripMargin

  implicit class OffsetDateTimeToSystemTimestamp(val offsetDateTime: OffsetDateTime) extends AnyVal {
    def toSystemTimestamp = Timestamp.valueOf(offsetDateTime.atZoneSameInstant(ZoneId.systemDefault).toLocalDateTime)
  }
}

trait MetadataStatement {
  def addKeyValue(key: String, value: Any): Unit
  def addEmptyValue(key: String): Unit
}

class MetadataStatementForWorkflow(preparedStatement: PreparedStatement, workflowId: String) extends MetadataStatement {
  import MetadataStatement._

  val offsetDateTimeFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ss.SSSZZZZZ")
  val logger = LoggerFactory.getLogger("LiquibaseMetadataMigration")
  val dawn = OffsetDateTime.of(0, 1, 1, 0, 0, 0, 0, ZoneOffset.UTC).toSystemTimestamp
  var batchSizeCounter: Int = 0

  private def metadataType(value: Any) = {
    value match {
      case WdlInteger(i) => "int"
      case WdlFloat(f) => "number"
      case WdlBoolean(b) => "boolean"
      case value: WdlValue => "string"
      case _: Int | Long => "int"
      case _: Double | Float => "number"
      case _: Boolean => "boolean"
      case _ =>"string"
    }
  }

  private def metadataValue(value: Any) = {
    value match {
      case v: WdlValue  => v.valueString
      case v => v.toString
    }
  }

  protected def setStatement() = {
    preparedStatement.setString(MetadataStatement.WorkflowIdIdx, workflowId)
    preparedStatement.setNull(MetadataStatement.CallFqnIdx, Types.VARCHAR)
    preparedStatement.setNull(MetadataStatement.CallIndexIdx, Types.INTEGER)
    preparedStatement.setNull(MetadataStatement.CallAttemptIdx, Types.INTEGER)
  }

  protected def addDataAndBatch(key: String, value: Any) = {
    preparedStatement.setString(MetadataStatement.KeyIdx, key)

    // Set the value and type
    value match {
      case null =>
        preparedStatement.setNull(MetadataStatement.ValueIdx, Types.VARCHAR)
        preparedStatement.setNull(MetadataStatement.ValueTypeIdx, Types.VARCHAR) // Null values have null type
      case _ =>
        preparedStatement.setString(MetadataStatement.ValueIdx, metadataValue(value))
        preparedStatement.setString(MetadataStatement.ValueTypeIdx, metadataType(value))
    }

    preparedStatement.addBatch()
  }

  private def add(key: String, value: Any, errorMessage: String) = try {
    setStatement()
    addDataAndBatch(key, value)
  } catch {
    case t: Throwable => logger.warn(errorMessage, t)
  }

  /** Adds a non-null value to the metadata journal. */
  override def addKeyValue(key: String, value: Any) = {
    if (value != null) {
      preparedStatement.setTimestamp(MetadataStatement.TimestampIdx, OffsetDateTime.now().toSystemTimestamp)
      add(key, value, s"Failed to migrate metadata value $value with key $key for workflow $workflowId")
    }
  }

  override def addEmptyValue(key: String): Unit = {
    preparedStatement.setTimestamp(MetadataStatement.TimestampIdx, dawn)
    add(key, null, s"Failed to add empty value with key $key for workflow $workflowId")
  }
}

class MetadataStatementForCall(preparedStatement: PreparedStatement, workflowId: String, callFqn: String, index: Int, attempt: Int) extends MetadataStatementForWorkflow(preparedStatement, workflowId) {
  override def setStatement() = {
    preparedStatement.setString(MetadataStatement.WorkflowIdIdx, workflowId)
    preparedStatement.setString(MetadataStatement.CallFqnIdx, callFqn)
    preparedStatement.setInt(MetadataStatement.CallIndexIdx, index)
    preparedStatement.setInt(MetadataStatement.CallAttemptIdx, attempt)
  }
}