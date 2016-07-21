package cromwell.database.slick

import java.io.{ByteArrayOutputStream, PrintStream}

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags._
import cromwell.database.DbmsTest
import cromwell.database.liquibase.DiffResultFilter
import liquibase.diff.output.DiffOutputControl
import liquibase.diff.output.changelog.DiffToChangeLog
import org.scalactic.StringNormalizations
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.{ExecutionContext, Future}
import scala.xml.{Elem, TopScope, XML}

class SlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures with StringNormalizations {

  behavior of "SlickDatabase"

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  it should "have the same liquibase and slick schema" in {
    for {
      liquibaseDatabase <- databaseForSchemaManager("liquibase").autoClosed
      slickDatabase <- databaseForSchemaManager("slick").autoClosed
    } {
      val diffResult = LiquibaseSchemaManager.compare(
        liquibaseDatabase.dataAccess.driver, liquibaseDatabase.database,
        slickDatabase.dataAccess.driver, slickDatabase.database)

      // TODO PBE get rid of this after the migration of #789 has run.
      val oldeTables = Seq(
        "EXECUTION_EVENT",
        "FAILURE_EVENT"
      )

      import DiffResultFilter._
      val diffFilters = StandardTypeFilters :+ UniqueIndexFilter
      val filteredDiffResult = diffResult
        .filterLiquibaseObjects
        .filterTableObjects(oldeTables)
        .filterChangedObjects(diffFilters)

      val totalChanged =
        filteredDiffResult.getChangedObjects.size +
        filteredDiffResult.getMissingObjects.size +
        filteredDiffResult.getUnexpectedObjects.size

      if (totalChanged > 0) {
        val outputStream = new ByteArrayOutputStream
        val printStream = new PrintStream(outputStream, true)
        val diffOutputControl = new DiffOutputControl(false, false, false, Array.empty)
        val diffToChangeLog = new DiffToChangeLog(filteredDiffResult, diffOutputControl)
        diffToChangeLog.print(printStream)
        val changeSetsScoped = XML.loadString(outputStream.toString) \ "changeSet" \ "_"
        val changeSets = changeSetsScoped map {
          case elem: Elem => elem.copy(scope = TopScope) // strip the namespace
          case other => other
        }
        fail(changeSets.mkString(
          "The following changes are in liquibase but not in slick:\n  ",
          "\n  ",
          "\nEither add the changes to slick or remove them from liquibase."))
      }
    }
  }


  it should "not deadlock" taggedAs PostMVP ignore {
    //    // Test based on https://github.com/kwark/slick-deadlock/blob/82525fc/src/main/scala/SlickDeadlock.scala
    //    val databaseConfig = ConfigFactory.parseString(
    //      s"""
    //         |db.url = "jdbc:hsqldb:mem:$${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
    //         |db.driver = "org.hsqldb.jdbcDriver"
    //         |db.numThreads = 2
    //         |driver = "slick.driver.HsqldbDriver$$"
    //         |""".stripMargin)
    //
    //    for {
    //      dataAccess <- (new SlickDatabase(databaseConfig) with DataAccess).autoClosed
    //    } {
    //      val futures = 1 to 20 map { _ =>
    //        val workflowId = WorkflowId.randomId()
    //        val workflowInfo = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = test1Sources)
    //        for {
    //          _ <- dataAccess.createWorkflow(workflowInfo, test1Sources, Nil, Nil, localBackend)
    //          _ <- dataAccess.getWorkflowExecutionAndAux(workflowInfo.id) map { result =>
    //            result.execution.workflowExecutionUuid should be(workflowId.toString)
    //          }
    //        } yield ()
    //      }
    //      Future.sequence(futures).futureValue(Timeout(10.seconds))
    //    }
  }

  def databaseForSchemaManager(schemaManager: String): SlickDatabase = {
    val databaseConfig = ConfigFactory.parseString(
      s"""
         |db.url = "jdbc:hsqldb:mem:$${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
         |db.driver = "org.hsqldb.jdbcDriver"
         |driver = "slick.driver.HsqldbDriver$$"
         |schema.manager = $schemaManager
         |""".stripMargin)
    new SlickDatabase(databaseConfig)
  }

  "SlickDatabase (main.hsqldb)" should behave like testWith("main.hsqldb")

  "SlickDatabase (test.mysql)" should behave like testWith("test.mysql")

  def testWith(configPath: String): Unit = {
    lazy val dataAccess: SlickDatabase = configPath match {
      case "main.hsqldb" => new SlickDatabase()
      case _ => new SlickDatabase(SlickDatabase.getDatabaseConfig(configPath))
    }

    it should "(if hsqldb) have transaction isolation mvcc" taggedAs DbmsTest in {
      import dataAccess.dataAccess.driver.api._

      val getProduct = SimpleDBIO[String](_.connection.getMetaData.getDatabaseProductName)
      //noinspection SqlDialectInspection
      val getHsqldbTx = sql"""SELECT PROPERTY_VALUE
                              FROM INFORMATION_SCHEMA.SYSTEM_PROPERTIES
                              WHERE PROPERTY_NAME = 'hsqldb.tx'""".as[String].head

      (for {
        product <- dataAccess.database.run(getProduct)
        _ <- product match {
          case "HSQL Database Engine" =>
            dataAccess.database.run(getHsqldbTx) map { hsqldbTx =>
              (hsqldbTx shouldEqual "mvcc") (after being lowerCased)
            }
          case _ => Future.successful(())
        }
      } yield ()).futureValue
    }
  }
}
