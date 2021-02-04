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
          containerOpt.foreach { _.stop }
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
          (metadataGenerator: MetadataGenerator) => {
            val metadata = metadataGenerator()
            (metadata.databaseMajorVersion, metadata.databaseMinorVersion) match {
              /*
              The docs say "If SLEEP() is interrupted, it returns 1."

              - https://mariadb.com/kb/en/sleep/

              Something changed in 10.3.26/10.4.16/10.5.7 that on a timeout started returning sqlState 70100 instead of
              a `1`. That triggers this line of code in the driver:
              - https://github.com/mariadb-corporation/mariadb-connector-j/blob/2.7.0/src/main/java/org/mariadb/jdbc/internal/util/exceptions/ExceptionFactory.java#L46-L48

              Skimming the change logs I can't quickly figure out what changed down in 10.3.26:
              - https://mariadb.com/kb/en/mariadb-10326-changelog/

              Note: 10.3.26/10.4.16/10.5.7 are all based on 10.2.35, but 10.2.35 is fine at the moment.

              If you find yourself here trying to return an IntErrorMessage(1) again for some new `latest`, you can also
              try reverting this commit to simplify the code. It's possible this is just a temporary bug introduced into
              10.3.X and the behavior will go back to prior behavior in 10.3.27 or some other future version.

              Also found this old discussion, but I wasn't able to fully understand it and determine what to expect:
              - https://jira.mariadb.org/browse/MDEV-10529?focusedCommentId=85479&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-85479

              Either way, we just need to make sure the sleep was interrupted and don't care if it was via an exception
              or returning `1`.
               */
              case (major, minor) if major >= 10 && minor >= 3 =>
                RegexErrorMessage(
                    """(\(conn=\d+\) )?Query execution was interrupted \(max_statement_time exceeded\)""".r
                )
              case _ => IntErrorMessage(1)
            }
          }
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
