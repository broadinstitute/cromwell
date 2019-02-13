package cromwell.database.migration.restart.table

import cromwell.core.simpleton.WomValueSimpleton._
import cromwell.database.migration.WdlTransformation._
import liquibase.database.jvm.JdbcConnection
import wdl.draft2.model.types.WdlFlavoredWomType
import wom.types.WomType
import wom.values.WomValue

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

    def buildJobStoreSimpletonEntries(name: String, womValue: WomValue, womType: WomType) = Option(womValue) match {
      case None => List(JobStoreSimpletonEntry(name, null, womType.stableName))
      case Some(_) => womValue.simplify(name) map { s =>
        JobStoreSimpletonEntry(s.simpletonKey, s.simpletonValue.valueString, s.simpletonValue.womType.stableName)
      }
    }

    while (results.next()) {
      val womType = WdlFlavoredWomType.fromDisplayString(results.getString(4))
      val rawValue = results.getString(3)
      val inflated = inflate(rawValue).get
      val womValue = coerceStringToWdl(inflated, womType)
      val name = results.getString(2)
      val jobStoreSimpletonEntry = buildJobStoreSimpletonEntries(name, womValue, womType)
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


