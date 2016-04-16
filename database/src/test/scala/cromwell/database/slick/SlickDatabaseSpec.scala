package cromwell.database.slick

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.database.liquibase.DiffResultFilter
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext

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
        slickDatabase.dataAccess.driver, slickDatabase.database,
        liquibaseDatabase.dataAccess.driver, liquibaseDatabase.database)

      import DiffResultFilter._
      val diffFilters = StandardTypeFilters :+ UniqueIndexFilter
      val filteredDiffResult = diffResult.filterLiquibaseObjects.filterChangedObjects(diffFilters)

      filteredDiffResult.getChangedObjects should be(empty)
      filteredDiffResult.getMissingObjects should be(empty)
      filteredDiffResult.getUnexpectedObjects should be(empty)
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
