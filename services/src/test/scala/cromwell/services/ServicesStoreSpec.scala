package cromwell.services

import java.io.{ByteArrayOutputStream, PrintStream}
import java.sql.Connection

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags._
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import liquibase.diff.DiffResult
import liquibase.diff.output.DiffOutputControl
import liquibase.diff.output.changelog.DiffToChangeLog
import org.hsqldb.persist.HsqlDatabaseProperties
import org.scalactic.StringNormalizations
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}
import slick.driver.JdbcProfile

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.xml._

class ServicesStoreSpec extends FlatSpec with Matchers with ScalaFutures with StringNormalizations {

  behavior of "ServicesStore"

  import ServicesStoreSpec._

  implicit val ec = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  it should "have the same liquibase and slick schema" in {
    for {
      liquibaseDatabase <- databaseForSchemaManager("liquibase").autoClosed
      slickDatabase <- databaseForSchemaManager("slick").autoClosed
    } {
      compare(
        liquibaseDatabase.dataAccess.driver, liquibaseDatabase.database,
        slickDatabase.dataAccess.driver, slickDatabase.database) { diffResult =>

        // TODO PBE get rid of this after the migration of #789 has run.
        val oldeTables = Seq(
          "EXECUTION",
          "EXECUTION_INFO",
          "EXECUTION_EVENT",
          "FAILURE_EVENT",
          "RUNTIME_ATTRIBUTES",
          "SYMBOL",
          "WORKFLOW_EXECUTION",
          "WORKFLOW_EXECUTION_AUX"
        )

        import cromwell.database.migration.liquibase.DiffResultFilter._
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
          val changeSets = changeSetsScoped map stripNodeScope
          fail(changeSets.mkString(
            "The following changes are in liquibase but not in slick:\n  ",
            "\n  ",
            "\nEither add the changes to slick or remove them from liquibase."))
        }
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

  "SlickDatabase (hsqldb)" should behave like testWith("database")

  "SlickDatabase (mysql)" should behave like testWith("database-test-mysql")

  def testWith(configPath: String): Unit = {
    import ServicesStore.EnhancedSqlDatabase

    lazy val databaseConfig = ConfigFactory.load.getConfig(configPath)
    lazy val dataAccess = new SlickDatabase(databaseConfig).initialized

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
          case HsqlDatabaseProperties.PRODUCT_NAME =>
            dataAccess.database.run(getHsqldbTx) map { hsqldbTx =>
              (hsqldbTx shouldEqual "mvcc") (after being lowerCased)
            }
          case _ => Future.successful(())
        }
      } yield ()).futureValue
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }
  }
}

object ServicesStoreSpec {
  // strip the namespace from elems and their children
  private def stripNodeScope(node: Node): Node = {
    node match {
      case elem: Elem => elem.copy(scope = TopScope, child = elem.child map stripNodeScope)
      case other => other
    }
  }

  private def databaseForSchemaManager(schemaManager: String): SlickDatabase = {
    val databaseConfig = ConfigFactory.parseString(
      s"""
         |db.url = "jdbc:hsqldb:mem:$${uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
         |db.driver = "org.hsqldb.jdbcDriver"
         |db.connectionTimeout = 3000
         |driver = "slick.driver.HsqldbDriver$$"
         |liquibase.updateSchema = false
         |""".stripMargin)
    val database = new SlickDatabase(databaseConfig)
    schemaManager match {
      case "liquibase" =>
        database withConnection LiquibaseUtils.updateSchema
      case "slick" =>
        SlickDatabase.createSchema(database)
    }
    database
  }

  private def compare[ReferenceProfile <: JdbcProfile, ComparisonProfile <: JdbcProfile, T]
  (referenceProfile: ReferenceProfile,
   referenceDatabase: ReferenceProfile#Backend#Database,
   comparisonProfile: ComparisonProfile,
   comparisonDatabase: ComparisonProfile#Backend#Database)(block: DiffResult => T)
  (implicit executor: ExecutionContext): T = {

    withConnections(referenceProfile, referenceDatabase, comparisonProfile, comparisonDatabase) {
      LiquibaseUtils.compare(_, _)(block)
    }
  }

  /**
    * Lends a connection to a block of code.
    *
    * @param profile  The slick jdbc profile for accessing the database.
    * @param database The database to use for the connection.
    * @param block    The block of code to run over the connection.
    * @tparam Profile The slick jdbc profile for accessing the database.
    * @tparam T       The return type of the block.
    * @return The return value of the block.
    */
  private def withConnection[Profile <: JdbcProfile, T](profile: Profile, database: Profile#Backend#Database)
                                                       (block: Connection => T): T = {
    /*
     TODO: Should this withConnection() method have a (implicit?) timeout parameter, that it passes on to Await.result?
     If we run completely asynchronously, nest calls to withConnection, and then call flatMap, the outer connection may
     already be closed before an inner block finishes running.
     */
    Await.result(database.run(profile.api.SimpleDBIO(context => block(context.connection))), Duration.Inf)
  }

  /**
    * Lends two connections to a block of code.
    *
    * @param profile1  The slick jdbc profile for accessing the first database.
    * @param database1 The database to use for the first connection.
    * @param profile2  The slick jdbc profile for accessing the second database.
    * @param database2 The database to use for the second connection.
    * @param block     The block of code to run over the first and second connections.
    * @tparam Profile1 The slick jdbc profile for accessing the first database.
    * @tparam Profile2 The slick jdbc profile for accessing the second database.
    * @tparam T        The return type of the block.
    * @return The return value of the block.
    */
  private def withConnections[Profile1 <: JdbcProfile, Profile2 <: JdbcProfile, T]
  (profile1: Profile1, database1: Profile1#Backend#Database, profile2: Profile2, database2: Profile2#Backend#Database)
  (block: (Connection, Connection) => T): T = {
    withConnection(profile1, database1) { connection1 =>
      withConnection(profile2, database2) { connection2 =>
        block(connection1, connection2)
      }
    }
  }
}
