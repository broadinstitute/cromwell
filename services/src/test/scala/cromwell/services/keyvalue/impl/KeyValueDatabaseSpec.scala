package cromwell.services.keyvalue.impl

import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.tables.JobKeyValueEntry
import cromwell.services.EngineServicesStore
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext

class KeyValueDatabaseSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {

  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  "SlickDatabase (hsqldb)" should behave like testWith("database")

  "SlickDatabase (mysql)" should behave like testWith("database-test-mysql")

  def testWith(configPath: String): Unit = {
    lazy val databaseConfig = ConfigFactory.load.getConfig(configPath)
    lazy val dataAccess = new EngineSlickDatabase(databaseConfig)
      .initialized(EngineServicesStore.EngineLiquibaseSettings)
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

    val keyValueEntryA2= JobKeyValueEntry(
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
  }
}
