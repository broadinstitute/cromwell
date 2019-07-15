package cromwell.services.database

import better.files._
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.Logger
import cromwell.core.WorkflowId
import cromwell.database.sql.joins.JobStoreJoin
import cromwell.database.sql.tables.JobStoreEntry
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}
import org.slf4j.LoggerFactory

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

/**
  * Tests for a regression in deadlocks within the Slick library.
  *
  * See cromwell.database.slick.SlickDatabase#actionThreadPool for more information.
  */
class SlickDeadlocksSpec extends FlatSpec with Matchers with ScalaFutures {

  implicit val executionContext = ExecutionContext.global

  behavior of "Slick Deadlocks"

  it should "not happen" in {
    val logger = Logger(LoggerFactory.getLogger(getClass.getName))
    // Test based on https://github.com/kwark/slick-deadlock/blob/82525fc/src/main/scala/SlickDeadlock.scala
    val databaseConfig = ConfigFactory.parseString(
      s"""|db.url = "jdbc:hsqldb:mem:$${uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
          |db.driver = "org.hsqldb.jdbcDriver"
          |db.connectionTimeout = 3000
          |db.numThreads = 2
          |profile = "slick.jdbc.HsqldbProfile$$"
          |""".stripMargin)
    logger.info("Initializing deadlock-test database")
    for {
      database <- DatabaseTestKit.initializedDatabaseFromConfig(EngineDatabaseType, databaseConfig).autoClosed
    } {
      logger.info(s"Initialized deadlock-test database: ${database.connectionDescription}")
      val futures = 1 to 20 map { count =>
        val workflowUuid = WorkflowId.randomId().toString
        val callFqn = "call.fqn"
        val jobIndex = 1
        val jobAttempt = 1
        val jobSuccessful = false
        val jobStoreEntry = JobStoreEntry(workflowUuid, callFqn, jobIndex, jobAttempt, jobSuccessful, None, None, None)
        val jobStoreJoins = Seq(JobStoreJoin(jobStoreEntry, Seq()))

        // NOTE: This test just needs to repeatedly read/write from a table that acts as a PK for a FK.
        logger.info(s"Creating deadlock-test future #$count with wf id $workflowUuid")
        val future = for {
          _ <- Future(logger.info(s"Starting build deadlock-test future #$count"))
          _ <- database.addJobStores(jobStoreJoins, 10)
          _ = logger.info(s"Added job stores for deadlock-test future #$count")
          queried <- database.queryJobStores(workflowUuid, callFqn, jobIndex, jobAttempt)
          _ = logger.info(s"Queried job stores for deadlock-test future #$count")
          _ = queried.get.jobStoreEntry.workflowExecutionUuid should be(workflowUuid)
          _ = logger.info(s"Successful check of deadlock-test future #$count")
        } yield ()
        logger.info(s"Created deadlock-test future #$count")
        future
      }

      logger.info("Futures created for deadlock-test. Waiting for sequence to complete…")
      Try(Future.sequence(futures).futureValue(Timeout(scaled(30.seconds)))) match {
        case Success(_) => logger.info("Passed deadlock-test")
        case Failure(throwable) =>
          logger.error("Failed deadlock-test")
          throw throwable
      }
    }
  }
}
