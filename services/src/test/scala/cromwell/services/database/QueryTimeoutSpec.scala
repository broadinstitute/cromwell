package cromwell.services.database

import cromwell.core.Tags.DbmsTest
import cromwell.database.slick.{MetadataSlickDatabase, SlickDatabase}
import org.scalactic.StringNormalizations
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers, PrivateMethodTester}

import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration._

class QueryTimeoutSpec extends FlatSpec with Matchers with ScalaFutures with StringNormalizations with PrivateMethodTester {

  // HSQL does not document a SLEEP() function, which is essential for this test
  // The functionality being tested is not relevant to an HSQL user, so the omission is probably acceptable
  val insomniacDatabases = Seq(HsqldbDatabaseSystem)

  // MariaDB driver does not appear to send the cancel command; there is not currently a requirement for this
  // feature to work for MariaDB customers
  val disobedientDatabases = Seq(MariadbDatabaseSystem)

  val databasesToTest = DatabaseSystem.All diff (insomniacDatabases union disobedientDatabases)

  val sleepCommands = Seq(
    "select sleep(10);",
    "select pg_sleep(10);"
  )

  val expectedErrors = Seq(
    "Statement cancelled due to timeout or client request",
    "ERROR: canceling statement due to user request"
  )

  for (((db, sleepCommand), errorMsg) <- databasesToTest zip sleepCommands zip expectedErrors) {
    behavior of s"${db.productName}"

    it should "fail with a timeout" taggedAs DbmsTest in {
      checkDatabaseSystem(db, sleepCommand, errorMsg)
    }
  }

  private def checkDatabaseSystem(databaseSystem: DatabaseSystem, sleepCommand: String, errorMsg: String) = {
    val database: SlickDatabase = DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)
    import database.dataAccess.driver.api._
    import slick.dbio.DBIO
    implicit val ec = ExecutionContext.global

    val testAdapter = new MetadataSlickDatabase(database.originalDatabaseConfig) {
      def runTestQuery[R](action: DBIO[R], timeout: Duration = Duration.Inf) = {
        runTransaction(action, timeout = timeout)
      }
    }

    import scala.language.reflectiveCalls // The compiler said to do this, I do not fully understand why

    val f: Future[Unit] = testAdapter.runTestQuery(sql"""#$sleepCommand""".as[Int], 5.seconds) map {
      _ =>
        fail("This query should have timed out!")
    } recover {
      case t =>
        assert(t.getMessage == errorMsg)
        ()
    }

    Await.result(f, 15.seconds)

    it should s"close the ${databaseSystem.productName} database" taggedAs DbmsTest in {
      liquibasedDatabase.close()
    }
  }
}
