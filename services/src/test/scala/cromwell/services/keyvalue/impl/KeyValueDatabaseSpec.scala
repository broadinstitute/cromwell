package cromwell.services.keyvalue.impl

import java.sql.SQLIntegrityConstraintViolationException

import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.sql.tables.JobKeyValueEntry
import cromwell.services.database._
import cromwell.services.keyvalue.impl.KeyValueDatabaseSpec._
import org.postgresql.util.PSQLException
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers, RecoverMethods}
import org.specs2.mock.Mockito

import scala.concurrent.{ExecutionContext, Future}

class KeyValueDatabaseSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito with RecoverMethods {

  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  DatabaseSystem.All foreach { databaseSystem =>
    behavior of s"KeyValueDatabase on ${databaseSystem.name}"

    lazy val dataAccess = DatabaseTestKit.initializedDatabaseFromSystem(EngineDatabaseType, databaseSystem)
    val workflowId = WorkflowId.randomId().toString
    val callFqn = "AwesomeWorkflow.GoodJob"

    val keyValueEntryA = JobKeyValueEntry(
      workflowExecutionUuid = workflowId,
      callFullyQualifiedName = callFqn,
      jobIndex = 0,
      jobAttempt = 1,
      storeKey = "myKeyA",
      storeValue = "myValueA"
    )

    val keyValueEntryA2 = JobKeyValueEntry(
      workflowExecutionUuid = workflowId,
      callFullyQualifiedName = callFqn,
      jobIndex = 0,
      jobAttempt = 1,
      storeKey = "myKeyA",
      storeValue = "myValueA2"
    )

    val keyValueEntryB = JobKeyValueEntry(
      workflowExecutionUuid = workflowId,
      callFullyQualifiedName = callFqn,
      jobIndex = 0,
      jobAttempt = 1,
      storeKey = "myKeyB",
      storeValue = "myValueB"
    )

    val wrongKeyValueEntryB = JobKeyValueEntry(
      workflowExecutionUuid = workflowId,
      callFullyQualifiedName = callFqn,
      jobIndex = 0,
      jobAttempt = 1,
      storeKey = "myKeyB",
      storeValue = null
    )

    it should "upsert and retrieve kv pairs correctly" taggedAs DbmsTest in {
      (for {
        // Just add A
        _ <- dataAccess.addJobKeyValueEntries(Seq(keyValueEntryA))
        valueA <- dataAccess.queryStoreValue(workflowId.toString, callFqn, 0, 1, "myKeyA")
        // Check that it's there
        _ = valueA shouldBe Some("myValueA")
        // Update A and add B in the same transaction
        _ <- dataAccess.addJobKeyValueEntries(Seq(keyValueEntryA2, keyValueEntryB))
        // A should have a new value
        valueA2 <- dataAccess.queryStoreValue(workflowId.toString, callFqn, 0, 1, "myKeyA")
        _ = valueA2 shouldBe Some("myValueA2")
        // B should also be there
        valueB <- dataAccess.queryStoreValue(workflowId.toString, callFqn, 0, 1, "myKeyB")
        _ = valueB shouldBe Some("myValueB")
      } yield ()).futureValue
    }

    it should "fail if one of the inserts fails" taggedAs DbmsTest in {
      val futureEx = recoverToExceptionIf[Exception] {
        dataAccess.addJobKeyValueEntries(Seq(keyValueEntryA, wrongKeyValueEntryB))
      }

      def verifyValues: Future[Unit] = for {
        // Values should not have changed
        valueA2 <- dataAccess.queryStoreValue(workflowId.toString, callFqn, 0, 1, "myKeyA")
        _ = valueA2 shouldBe Some("myValueA2")
        // B should also be there
        valueB <- dataAccess.queryStoreValue(workflowId.toString, callFqn, 0, 1, "myKeyB")
        _ = valueB shouldBe Some("myValueB")
      } yield ()

      (futureEx map { ex =>
        ex.getMessage should fullyMatch regex getFailureRegex(databaseSystem)
        ex.getClass should be(getFailureClass(databaseSystem))
      }).flatMap(_ => verifyValues).futureValue
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }
  }
}

object KeyValueDatabaseSpec {
  private def getFailureRegex(databaseSystem: DatabaseSystem): String = {
    databaseSystem.platform match {
      case HsqldbDatabasePlatform =>
        "integrity constraint violation: NOT NULL check constraint; SYS_CT_10591 table: JOB_KEY_VALUE_ENTRY column: STORE_VALUE"
      case MariadbDatabasePlatform => """\(conn=\d+\) Column 'STORE_VALUE' cannot be null"""
      case MysqlDatabasePlatform => "Column 'STORE_VALUE' cannot be null"
      case PostgresqlDatabasePlatform => """ERROR: null value in column "STORE_VALUE" violates not-null constraint"""
    }
  }

  private def getFailureClass(databaseSystem: DatabaseSystem): Class[_ <: Exception] = {
    databaseSystem.platform match {
      case HsqldbDatabasePlatform => classOf[SQLIntegrityConstraintViolationException]
      case MariadbDatabasePlatform => classOf[SQLIntegrityConstraintViolationException]
      case MysqlDatabasePlatform => classOf[SQLIntegrityConstraintViolationException]
      case PostgresqlDatabasePlatform => classOf[PSQLException]
    }
  }
}
