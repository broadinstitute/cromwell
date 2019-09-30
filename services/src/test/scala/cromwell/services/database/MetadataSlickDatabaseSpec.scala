package cromwell.services.database

import java.time.OffsetDateTime

import cromwell.core.Tags.DbmsTest
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.concurrent.ExecutionContext

class MetadataSlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures {

  DatabaseSystem.All foreach { databaseSystem =>

    implicit val ec = ExecutionContext.global

    behavior of s"MetadataSlickDatabase on ${databaseSystem.name}"

    lazy val database: MetadataSlickDatabase with TestSlickDatabase = DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)
    import database.dataAccess.driver.api._

    import cromwell.database.migration.metadata.table.symbol.MetadataStatement.OffsetDateTimeToSystemTimestamp
    val now = OffsetDateTime.now().toSystemTimestamp

    it should "set up the test data" taggedAs DbmsTest in {
      database.runTestTransaction(
        database.dataAccess.metadataEntries ++= Seq(
          MetadataEntry("workflow id: 3 to delete, 1 label", None, None, None, "someKey", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("workflow id: 3 to delete, 1 label", None, None, None, "someKey", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("workflow id: 3 to delete, 1 label", None, None, None, "someKey", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("workflow id: 3 to delete, 1 label", None, None, None, "labels:dontDeleteMe", None, None, OffsetDateTime.now().toSystemTimestamp, None),

          MetadataEntry("workflow id: I am a root workflow with a subworkflow", None, None, None, "labels:dontDeleteMe", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("workflow id: I am a root workflow with a subworkflow", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("workflow id: I am the subworkflow", None, None, None, "labels:dontDeleteMe", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("workflow id: I am the subworkflow", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),

          MetadataEntry("nested subworkflows: root", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("nested subworkflows: first nesting", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("nested subworkflows: second nesting", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),
          MetadataEntry("nested subworkflows: third nesting nesting", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),
        )
      ).futureValue(Timeout(10.seconds))

      database.runTestTransaction(
        database.dataAccess.workflowMetadataSummaryEntries ++= Seq(
          WorkflowMetadataSummaryEntry("workflow id: I am not a root workflow", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("I have a parent"), Option("I have a parent"), None),
          WorkflowMetadataSummaryEntry("workflow id: 3 to delete, 1 label", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), None, None, None),

          WorkflowMetadataSummaryEntry("workflow id: I am a root workflow with a subworkflow", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), None, None, None),
          WorkflowMetadataSummaryEntry("workflow id: I am the subworkflow", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("workflow id: I am a root workflow with a subworkflow"), Option("workflow id: I am a root workflow with a subworkflow"), None),

          WorkflowMetadataSummaryEntry("workflow id: I am still running!", Option("workflow name"), Option("Running"), Option(now), Option(now), Option(now), Option("workflow id: I am a root workflow with a subworkflow"), Option("workflow id: I am a root workflow with a subworkflow"), None),
          WorkflowMetadataSummaryEntry("workflow id: I inexplicably do not have a workflow status", Option("workflow name"), None, Option(now), Option(now), Option(now), Option("workflow id: I am a root workflow with a subworkflow"), Option("workflow id: I am a root workflow with a subworkflow"), None),

          WorkflowMetadataSummaryEntry("nested subworkflows: root", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), None, None, None),
          WorkflowMetadataSummaryEntry("nested subworkflows: first nesting", Option("workflow name"), None, Option(now), Option(now), Option(now), Option("nested subworkflows: root"), Option("nested subworkflows: root"), None),
          WorkflowMetadataSummaryEntry("nested subworkflows: second nesting", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("nested subworkflows: first nesting"), Option("nested subworkflows: root"), None),
          WorkflowMetadataSummaryEntry("nested subworkflows: third nesting nesting", Option("workflow name"), None, Option(now), Option(now), Option(now), Option("nested subworkflows: third nesting nesting"), Option("nested subworkflows: root"), None),
        )
      ).futureValue(Timeout(10.seconds))
    }

    it should "error when deleting a root workflow that does not exist" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("does not exist")

      delete.failed.futureValue(Timeout(10.seconds)).getMessage should be("""[Carbonite metadata deletion] Failed with non-existent root workflow "does not exist"""")
    }

    it should "error when calling delete on a subworkflow ID" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("workflow id: I am not a root workflow")
      delete.failed.futureValue(Timeout(10.seconds)).getMessage should be("""[Carbonite metadata deletion] Failed because workflow is not root: "workflow id: I am not a root workflow"""")
    }

    it should "error when trying to delete a workflow that has subworkflows itself" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("nested subworkflows: second nesting")
      delete.failed.futureValue(Timeout(10.seconds)).getMessage should be("""[Carbonite metadata deletion] Failed because workflow is not root: "nested subworkflows: second nesting"""")
    }

    it should "error when trying to delete a still running workflow" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("workflow id: I am still running!")
      delete.failed.futureValue(Timeout(10.seconds)).getMessage should be("""[Carbonite metadata deletion] Failed with non-terminal summary status "Running" for workflow "workflow id: I am still running!"""")
    }

    it should "error when trying to delete a workflow with missing workflow status" taggedAs DbmsTest in {
      // I don't know how this would happen but it is technically nullable, and pretty important to not screw up
      val delete = database.deleteNonLabelMetadataForWorkflow("workflow id: I inexplicably do not have a workflow status")
      delete.failed.futureValue(Timeout(10.seconds)).getMessage should be("""[Carbonite metadata deletion] Failed because summary status unexpectedly empty for workflow "workflow id: I inexplicably do not have a workflow status"""")
    }

    it should "delete the right number of rows for a root workflow without subworkflows" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("workflow id: 3 to delete, 1 label")
      delete.futureValue(Timeout(10.seconds)) should be(3)
    }

    it should "delete the right number of rows for a root workflow with subworkflows" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("workflow id: I am a root workflow with a subworkflow")
      delete.futureValue(Timeout(10.seconds)) should be(2)
    }

    it should "delete the right number of rows for a nested subworkflow" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflow("nested subworkflows: root")
      delete.futureValue(Timeout(10.seconds)) should be(4)
    }

    it should "clean up & close the database" taggedAs DbmsTest in {
      // Not relevant in Travis where all state gets nuked but useful for testing locally
      database.runTestTransaction(database.dataAccess.metadataEntries.delete).futureValue(Timeout(10.seconds))
      database.runTestTransaction(database.dataAccess.workflowMetadataSummaryEntries.delete).futureValue(Timeout(10.seconds))

      database.close()
    }

  }


}
