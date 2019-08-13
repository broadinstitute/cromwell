package cromwell.services.database

import better.files._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.Tags.DbmsTest
import cromwell.database.slick.SlickDatabase
import cromwell.database.slick.tables.DataAccessComponent
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.Future
import scala.concurrent.duration._

class QueryTimeoutSpec extends FlatSpec with Matchers with ScalaFutures {

  // HSQL does not document a SLEEP() function, which is essential for this test
  // The functionality being tested is not relevant to an HSQL user, so the omission is probably acceptable
  val insomniacDatabases = Seq(HsqldbDatabaseSystem)

  val databasesToTest = DatabaseSystem.All diff insomniacDatabases

  val sleepCommands = Seq(
    "select sleep(10);",
    "select sleep(10);",
    "select pg_sleep(10);"
  )

  val expectedErrors = Seq(
    Right(Option(1)),
    Left("Statement cancelled due to timeout or client request"),
    Left("ERROR: canceling statement due to user request"),
  )

  for (((db, sleepCommand), errorMsg) <- databasesToTest zip sleepCommands zip expectedErrors) {
    behavior of s"${db.productName}"

    it should "fail with a timeout" taggedAs DbmsTest in {
      checkDatabaseSystem(db, sleepCommand, errorMsg)
    }
  }

  private def checkDatabaseSystem(databaseSystem: DatabaseSystem,
                                  sleepCommand: String,
                                  errorEither: Either[String, Option[Int]]): Unit = {

    for {
      testDatabase <- new TestDatabase(ConfigFactory.load.getConfig(databaseSystem.configPath)).autoClosed
    } {
      import testDatabase.dataAccess.driver.api._

      //noinspection SqlDialectInspection
      val future = testDatabase.runTestTransaction(sql"""#$sleepCommand""".as[Int].headOption, 5.seconds)
      errorEither match {
        case Left(message) => future.failed.futureValue(Timeout(10.seconds)).getMessage should be(message)
        case Right(optionResult) => future.futureValue(Timeout(10.seconds)) should be(optionResult)
      }
    }
  }
}

class TestDatabase(config: Config) extends SlickDatabase(config) {
  override lazy val dataAccess: DataAccessComponent = new DataAccessComponent {
    override lazy val driver = slickConfig.profile
    override lazy val schema = driver.DDL("", "")
  }

  import dataAccess.driver.api._

  /**
    * Run a RepeatableRead transaction directly on the database with an optional timeout.
    */
  def runTestTransaction[R](action: DBIO[R], timeout: Duration = Duration.Inf): Future[R] = {
    runTransaction(action, timeout = timeout)
  }
}
