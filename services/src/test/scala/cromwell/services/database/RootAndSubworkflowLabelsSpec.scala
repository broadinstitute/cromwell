package cromwell.services.database

import java.time.OffsetDateTime

import cromwell.core.Tags.DbmsTest
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.database.sql.tables.{CustomLabelEntry, WorkflowMetadataSummaryEntry}
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class RootAndSubworkflowLabelsSpec extends FlatSpec with Matchers with ScalaFutures {
  implicit val ec = ExecutionContext.global

  DatabaseSystem.All foreach { databaseSystem =>
    behavior of s"MetadataSlickDatabase on ${databaseSystem.name}"

    lazy val database: MetadataSlickDatabase with TestSlickDatabase = DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)
    import database.dataAccess.driver.api._

    import cromwell.database.migration.metadata.table.symbol.MetadataStatement.OffsetDateTimeToSystemTimestamp
    val now = OffsetDateTime.now().toSystemTimestamp
    println(now)

    it should "set up the test data" taggedAs DbmsTest in {
      database.runTestTransaction(
        {
          val root = summary(name = "root")
          val branch = summary(name = "branch", parent = Option("root"), root = Option("root"))
          val leaf = summary(name = "leaf", parent = Option("branch"), root = Option("root"))
          val random = summary(name = "random")
          database.dataAccess.workflowMetadataSummaryEntries ++= Seq(root, branch, leaf, random)
        } andThen {
          val root = label("root")
          // intentionally not labeling the branch as a negative test
          val leaf = label("leaf")
          val random = label("random")
          database.dataAccess.customLabelEntries ++= Seq(root, /* branch, */ leaf, random)
        }
      ).futureValue(Timeout(10.seconds))
    }

    def summary(name: String, parent: Option[String] = None, root: Option[String] = None): WorkflowMetadataSummaryEntry = {
      WorkflowMetadataSummaryEntry(
        workflowExecutionUuid = name,
        workflowStatus = Option("Succeeded"),
        workflowName = Option(name),
        startTimestamp = None,
        endTimestamp = None,
        submissionTimestamp = None,
        parentWorkflowExecutionUuid = parent,
        rootWorkflowExecutionUuid = root,
        metadataArchiveStatus = None
      )
    }

    def label(uuid: String): CustomLabelEntry = {
      CustomLabelEntry(customLabelKey = "key", customLabelValue = uuid, workflowExecutionUuid = uuid)
    }

    it should "query root and subworkflow labels correctly" taggedAs DbmsTest in {
      database.getRootAndSubworkflowLabels("root").
        futureValue(Timeout(10.seconds)) shouldBe Map(
        "root" -> Map("key" -> "root"),
        "leaf" -> Map("key" -> "leaf")
      )
    }
  }
}
