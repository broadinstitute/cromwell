package cromwell.database.migration.restart.table

import java.io.{ByteArrayInputStream, IOException}
import java.util.zip.GZIPInputStream

import cromwell.core.simpleton.WdlValueSimpleton._
import liquibase.change.custom.CustomTaskChange
import liquibase.database.Database
import liquibase.database.jvm.JdbcConnection
import liquibase.exception.ValidationErrors
import liquibase.resource.ResourceAccessor
import org.apache.commons.codec.binary.Base64
import org.apache.commons.io.IOUtils
import wdl4s.types.{WdlPrimitiveType, WdlType}

import scala.language.postfixOps
import scala.util.Try


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

  private def inflate(value: String) = {
    Try {
      IOUtils.toString(new GZIPInputStream(new ByteArrayInputStream(Base64.decodeBase64(value))))
    } recover {
      case e: IOException => value
    } get
  }

  protected def doMigration(connection: JdbcConnection): Unit = {
    val query = connection.createStatement()
    lazy val insert = connection.prepareStatement(InsertJobStoreSimpleton)
    query.execute(QueryOutputsForDoneCallsInRunningWorkflows)
    val results = query.getResultSet
    while (results.next()) {
      val wdlType = WdlType.fromWdlString(results.getString(4))
      val rawValue = results.getString(3)
      val inflated = inflate(rawValue)
      val wdlValue = wdlType match {
        case p: WdlPrimitiveType => p.coerceRawValue(inflated).get
        case o => o.fromWdlString(inflated)
      }
      val name = results.getString(2)
      val simpletons = wdlValue.simplify(name)
      simpletons foreach { s =>
        insert.setInt(1, results.getInt(1))
        insert.setString(2, s.simpletonKey)
        insert.setString(3, s.simpletonValue.toWdlString)
        insert.setString(4, s.simpletonValue.wdlType.toWdlString)
        insert.execute()
      }
    }
  }
}


