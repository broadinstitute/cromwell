package cromwell.database.migration.restart.table

import cromwell.core.simpleton.WdlValueSimpleton._
import cromwell.database.migration.WdlTransformation._
import liquibase.database.jvm.JdbcConnection
import wdl4s.types.WdlType
import wdl4s.values.WdlValue


class JobStoreSimpletonMigration extends AbstractRestartMigration {
  override val description = "WORKFLOW_EXECUTION + EXECUTION + SYMBOL + JOB_STORE -> JOB_STORE_RESULT_SIMPLETON"

  // GOTC (substituting COUNT(*) for the projection): 1 row in set (5.22 sec)
  private val QueryOutputsForDoneCallsInRunningWorkflows =
  """
      SELECT
        js.JOB_STORE_ID,   -- 1
        s.NAME,            -- 2
        s.WDL_VALUE,       -- 3
        s.WDL_TYPE         -- 4
      FROM EXECUTION e
        JOIN WORKFLOW_EXECUTION we
          ON we.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID
        JOIN SYMBOL s
          ON s.WORKFLOW_EXECUTION_ID = e.WORKFLOW_EXECUTION_ID
        JOIN JOB_STORE js
          ON js.WORKFLOW_UUID     = we.WORKFLOW_EXECUTION_UUID AND
             js.CALL_FQN          = e.CALL_FQN AND
             js.JOB_SCATTER_INDEX = e.IDX AND
             js.JOB_RETRY_ATTEMPT = e.ATTEMPT
      WHERE
        s.IO = 'OUTPUT' AND
        s.SCOPE = js.CALL_FQN AND
        s.INDEX = js.JOB_SCATTER_INDEX
  """

  private val InsertJobStoreSimpleton =
    """
      INSERT INTO JOB_STORE_RESULT_SIMPLETON(
        JOB_STORE_ID, SIMPLETON_KEY, SIMPLETON_VALUE, WDL_TYPE)
          VALUES(?, ?, ?, ?)
    """

  protected def doMigration(connection: JdbcConnection): Unit = {
    val query = connection.createStatement()
    lazy val insert = connection.prepareStatement(InsertJobStoreSimpleton)
    query.execute(QueryOutputsForDoneCallsInRunningWorkflows)
    val results = query.getResultSet

    case class JobStoreSimpletonEntry(name: String, valueString: String, typeString: String)

    def buildJobStoreSimpletonEntries(name: String, wdlValue: WdlValue, wdlType: WdlType) = Option(wdlValue) match {
      case None => List(JobStoreSimpletonEntry(name, null, wdlType.toWdlString))
      case Some(v) => wdlValue.simplify(name) map { s =>
        JobStoreSimpletonEntry(s.simpletonKey, s.simpletonValue.toWdlString, s.simpletonValue.wdlType.toWdlString)
      }
    }

    while (results.next()) {
      val wdlType = WdlType.fromWdlString(results.getString(4))
      val rawValue = results.getString(3)
      val inflated = inflate(rawValue).get
      val wdlValue = coerceStringToWdl(inflated, wdlType)
      val name = results.getString(2)
      val jobStoreSimpletonEntry = buildJobStoreSimpletonEntries(name, wdlValue, wdlType)
      jobStoreSimpletonEntry foreach { e =>
        insert.setInt(1, results.getInt(1))
        insert.setString(2, e.name)
        insert.setString(3, e.valueString)
        insert.setString(4, e.typeString)
        insert.execute()
      }
    }
  }
}


