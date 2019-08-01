package cromwell.services.database

import cromwell.core.Tags.DbmsTest
import cromwell.database.slick.{MetadataSlickDatabase, SlickDatabase}
import org.scalactic.StringNormalizations
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers, PrivateMethodTester}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Success}

class QueryTimeoutSpec extends FlatSpec with Matchers with ScalaFutures with StringNormalizations with PrivateMethodTester {

  // HSQL does not document a SLEEP() function, which is essential for this test
  // The functionality being tested is not relevant to an HSQL user, so the omission is probably acceptable
  val insomniacDatabases = Seq(HsqldbDatabaseSystem)

  val databasesToTest = DatabaseSystem.All diff insomniacDatabases

  val sleepCommands = Seq(
    "select sleep(10);",
    "select sleep(10);",
    "select pg_sleep(10);"
  )

  for ((db, sleepCommand) <- databasesToTest zip sleepCommands) {
    behavior of s"${db.productName}"

    it should "fail with a timeout" taggedAs DbmsTest in {
      checkDatabaseSystem(db, sleepCommand)
    }
  }

  private def checkDatabaseSystem(databaseSystem: DatabaseSystem, sleepCommand: String) = {
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

    testAdapter.runTestQuery(sqlu"$sleepCommand", 5.seconds) onComplete {
      case Success(a) =>
        println(a)
        fail // TODO: make sure we're failing for the right reason, i.e. timeout, not that the DB is unavailable
      case Failure(t: Throwable) =>
        println(t)
        succeed
    }
  }
}
