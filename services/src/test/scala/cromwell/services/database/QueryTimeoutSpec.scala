package cromwell.services.database

import better.files._
import com.dimafeng.testcontainers.Container
import common.assertion.CromwellTimeoutSpec
import cromwell.core.Tags.DbmsTest
import cromwell.services.database.QueryTimeoutSpec._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.util.matching.Regex

/**
  * Tests that when we pass a timeout into databases that it actually times out.
  *
  * Runs a DBMS specific `SLEEP(10)` and tells the database driver to timeout the query after 5 seconds.
  * Does not test database systems that don't support some sort of `SLEEP()` SQL function.
  */
class QueryTimeoutSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with ScalaFutures {

  DatabaseSystem.All foreach { databaseSystem =>
    testOption(databaseSystem) foreach {
      case (sleepCommand, errorMessageGenerator) =>
        behavior of s"Query timeouts on ${databaseSystem.name}"

        val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

        it should "start container if required" taggedAs DbmsTest in {
          containerOpt.foreach { _.start }
        }

        it should "fail with a timeout" taggedAs DbmsTest in {
          checkDatabaseSystem(containerOpt, databaseSystem, sleepCommand, errorMessageGenerator)
        }

        it should "stop container if required" taggedAs DbmsTest in {
          containerOpt.foreach { _.stop() }
        }
    }
  }

  private def checkDatabaseSystem(containerOpt: Option[Container],
                                  databaseSystem: DatabaseSystem,
                                  sleepCommand: String,
                                  errorMessageGenerator: ErrorMessageGenerator): Unit = {
    for {
      testDatabase <- DatabaseTestKit.schemalessDatabaseFromContainerOptAndSystem(containerOpt, databaseSystem).autoClosed
    } {
      import testDatabase.dataAccess.driver.api._

      def metadataGenerator(): ConnectionMetadata = DatabaseTestKit.connectionMetadata(testDatabase)

      val errorMessage = errorMessageGenerator(metadataGenerator)

      //noinspection SqlDialectInspection
      val future = testDatabase.runTestTransaction(sql"""#$sleepCommand""".as[Int].headOption, timeout = 5.seconds)
      errorMessage match {
        case IntErrorMessage(result) => future.futureValue(Timeout(10.seconds)) should be(Option(result))
        case StringErrorMessage(message) => future.failed.futureValue(Timeout(10.seconds)).getMessage should be(message)
        case RegexErrorMessage(pattern) =>
          future.failed.futureValue(Timeout(10.seconds)).getMessage should fullyMatch regex pattern
      }
    }
  }

  private def testOption(databaseSystem: DatabaseSystem): Option[(String, ErrorMessageGenerator)] = {
    databaseSystem.platform match {
      case HsqldbDatabasePlatform =>
        // HSQL does not document a SLEEP() function, which is essential for this test
        // The functionality being tested is not relevant to an HSQL user, so the omission is probably acceptable
        None
      case MariadbDatabasePlatform =>
        Option((
          "select sleep(10);",
          _ => RegexErrorMessage("""(\(conn=\d+\) )?Query execution was interrupted \(max_statement_time exceeded\)""".r)
        ))
      case MysqlDatabasePlatform =>
        Option((
          "select sleep(10);",
          _ => StringErrorMessage("Statement cancelled due to timeout or client request"),
        ))
      case PostgresqlDatabasePlatform =>
        Option((
          "select pg_sleep(10);",
          _ => StringErrorMessage("ERROR: canceling statement due to user request"),
        ))
    }
  }
}

object QueryTimeoutSpec {
  type MetadataGenerator = () => ConnectionMetadata
  type ErrorMessageGenerator = MetadataGenerator => ErrorMessage

  sealed trait ErrorMessage

  case class IntErrorMessage(value: Int) extends ErrorMessage

  case class StringErrorMessage(value: String) extends ErrorMessage

  case class RegexErrorMessage(value: Regex) extends ErrorMessage
}
