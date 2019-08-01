package cromwell.services.database

import cromwell.core.Tags.DbmsTest
import cromwell.database.slick.MetadataSlickDatabase
import org.scalactic.StringNormalizations
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers, PrivateMethodTester}
import slick.jdbc.TransactionIsolation

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.util.{Failure, Success}

class QueryTimeoutSpec extends FlatSpec with Matchers with ScalaFutures with StringNormalizations with PrivateMethodTester {

  // HSQL does not document a SLEEP() function, which is essential for this test
  // The functionality being tested is not relevant to an HSQL user, so the omission is probably acceptable
  val insomniacDatabases = Seq(HsqldbDatabaseSystem)

  val databasesToTest = DatabaseSystem.All diff insomniacDatabases

  val sleepCommands = Seq(
    "select pg_sleep(10);",
    "select sleep(10);",
    "select sleep(10);"
  )

  for ((db, sleepCommand) <- databasesToTest zip sleepCommands) {
    behavior of s"${db.productName}"

    it should "fail with a timeout" taggedAs DbmsTest in {
      checkDatabaseSystem(PostgresqlDatabaseSystem, sleepCommand)
    }
  }

  private def checkDatabaseSystem(databaseSystem: DatabaseSystem, sleepCommand: String) = {
    val database: MetadataSlickDatabase = DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)
    import database.dataAccess.driver.api._

    val runTransactionMethod = PrivateMethod[MetadataSlickDatabase]('runTransaction)

    val transactionFuture = runTransactionMethod(sqlu"$sleepCommand", TransactionIsolation.RepeatableRead, 5.seconds).asInstanceOf[Future[_]]
    implicit val ec = ExecutionContext.global

    transactionFuture onComplete {
      case Success(_) =>
        fail // TODO: make sure we're failing for the right reason, i.e. timeout, not that the DB is unavailable
      case Failure(_) =>
        succeed
    }
  }
}
