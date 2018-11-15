package cromwell.engine.workflow.workflowstore

import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.database.slick.EngineSlickDatabase
import cromwell.services.EngineServicesStore
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class WorkflowStoreSlickSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(10, Seconds)), scaled(Span(100, Millis)))

  "SlickDatabase (hsqldb)" should behave like testWith("database")

  "SlickDatabase (mysql)" should behave like testWith("database-test-mysql")

  def testWith(configPath: String): Unit = {
    lazy val databaseConfig = ConfigFactory.load.getConfig(configPath)
    val sourceFilesCollection = NonEmptyList.of(WorkflowSourceFilesCollection(Option("sample"), None, None, None, None, "input", "option", "string", None, workflowOnHold = true, Seq.empty))

    it should "honor the onHold flag" taggedAs DbmsTest in {
      val workflowStore: WorkflowStore = SqlWorkflowStore(new EngineSlickDatabase(databaseConfig)
        .initialized(EngineServicesStore.EngineLiquibaseSettings))

      (for {
        submissionResponses <- workflowStore.add(sourceFilesCollection)
        workflowCount <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second)
        _ = workflowCount.size shouldBe 0
        _ <- workflowStore.switchOnHoldToSubmitted(submissionResponses.head.id)
        workflowCount2 <- workflowStore.fetchStartableWorkflows(10, "A00", 1.second)
        _ = workflowCount2.size shouldBe 1
      } yield ()).futureValue
    }
  }
}
