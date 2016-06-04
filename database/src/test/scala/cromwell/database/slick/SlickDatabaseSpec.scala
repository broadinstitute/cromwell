package cromwell.database.slick

import java.io.{ByteArrayOutputStream, PrintStream}

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.database.liquibase.DiffResultFilter
import liquibase.diff.output.DiffOutputControl
import liquibase.diff.output.changelog.DiffToChangeLog
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext
import scala.xml.{Elem, TopScope, XML}

class SlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures {

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
}
