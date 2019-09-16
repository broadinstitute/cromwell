package cromwell.services.database

import better.files._
import cromwell.core.Tags.DbmsTest
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._

class QueryTimeoutSpec extends FlatSpec with Matchers with ScalaFutures {

  DatabaseSystem.All foreach { databaseSystem =>
    testOption(databaseSystem) foreach {
      case (sleepCommand, errorMessage) =>
        behavior of s"Query timeouts on ${databaseSystem.name}"

        it should "fail with a timeout" taggedAs DbmsTest in {
          checkDatabaseSystem(databaseSystem, sleepCommand, errorMessage)
        }
    }
  }

  private def checkDatabaseSystem(databaseSystem: DatabaseSystem,
                                  sleepCommand: String,
                                  errorEither: Either[String, Option[Int]]): Unit = {
    for {
      testDatabase <- DatabaseTestKit.schemalessDatabaseFromSystem(databaseSystem).autoClosed
    } {
      import testDatabase.dataAccess.driver.api._

      //noinspection SqlDialectInspection
      val future = testDatabase.runTestTransaction(sql"""#$sleepCommand""".as[Int].headOption, timeout = 5.seconds)
      errorEither match {
        case Left(message) => future.failed.futureValue(Timeout(10.seconds)).getMessage should be(message)
        case Right(optionResult) => future.futureValue(Timeout(10.seconds)) should be(optionResult)
      }
    }
  }

  private def testOption(databaseSystem: DatabaseSystem): Option[(String, Either[String, Option[Int]])] = {
    databaseSystem.platform match {
      case HsqldbDatabasePlatform =>
        // HSQL does not document a SLEEP() function, which is essential for this test
        // The functionality being tested is not relevant to an HSQL user, so the omission is probably acceptable
        None
      case MariadbDatabasePlatform =>
        Option((
          "select sleep(10);",
          Right(Option(1)),
        ))
      case MysqlDatabasePlatform =>
        Option((
          "select sleep(10);",
          Left("Statement cancelled due to timeout or client request"),
        ))
      case PostgresqlDatabasePlatform =>
        Option((
          "select pg_sleep(10);",
          Left("ERROR: canceling statement due to user request"),
        ))
    }
  }
}
