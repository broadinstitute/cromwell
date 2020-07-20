package cromwell.services.database

import java.time.OffsetDateTime

import com.dimafeng.testcontainers.Container
import cromwell.core.Tags.DbmsTest
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.database.migration.metadata.table.symbol.MetadataStatement._
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.database.slick.MetadataSlickDatabase.SummarizationPartitionedMetadata
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.services.metadata.CallMetadataKeys
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class MetadataSlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures {

  DatabaseSystem.All foreach { databaseSystem =>

    implicit val ec = ExecutionContext.global

    behavior of s"MetadataSlickDatabase on ${databaseSystem.name}"

    val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

    lazy val database = DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, MetadataDatabaseType, databaseSystem)
    import cromwell.database.migration.metadata.table.symbol.MetadataStatement.OffsetDateTimeToSystemTimestamp
    import database.dataAccess.driver.api._
    val now = OffsetDateTime.now().toSystemTimestamp

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach { _.start }
    }

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
          MetadataEntry("nested subworkflows: third nesting", None, None, None, "please do delete me", None, None, OffsetDateTime.now().toSystemTimestamp, None),
        )
      ).futureValue(Timeout(10.seconds))

      database.runTestTransaction(
        database.dataAccess.workflowMetadataSummaryEntries ++= Seq(
          WorkflowMetadataSummaryEntry("workflow id: 3 to delete, 1 label", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), None, None, None, None),

          WorkflowMetadataSummaryEntry("workflow id: I am a root workflow with a subworkflow", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), None, None, None, None),
          WorkflowMetadataSummaryEntry("workflow id: I am the subworkflow", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("workflow id: I am a root workflow with a subworkflow"), Option("workflow id: I am a root workflow with a subworkflow"), None, None),

          WorkflowMetadataSummaryEntry("nested subworkflows: root", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), None, None, None, None),
          WorkflowMetadataSummaryEntry("nested subworkflows: first nesting", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("nested subworkflows: root"), Option("nested subworkflows: root"), None, None),
          WorkflowMetadataSummaryEntry("nested subworkflows: second nesting", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("nested subworkflows: first nesting"), Option("nested subworkflows: root"), None, None),
          WorkflowMetadataSummaryEntry("nested subworkflows: third nesting", Option("workflow name"), Option("Succeeded"), Option(now), Option(now), Option(now), Option("nested subworkflows: second nesting"), Option("nested subworkflows: root"), None, None),
        )
      ).futureValue(Timeout(10.seconds))
    }

    it should "delete the right number of rows for a root workflow without subworkflows" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflowAndUpdateArchiveStatus("workflow id: 3 to delete, 1 label", None)
      delete.futureValue(Timeout(10.seconds)) should be(3)
    }

    it should "delete the right number of rows for a root workflow with subworkflows" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflowAndUpdateArchiveStatus("workflow id: I am a root workflow with a subworkflow", None)
      delete.futureValue(Timeout(10.seconds)) should be(2)
    }

    it should "delete the right number of rows for a nested subworkflow" taggedAs DbmsTest in {
      val delete = database.deleteNonLabelMetadataForWorkflowAndUpdateArchiveStatus("nested subworkflows: root", None)
      delete.futureValue(Timeout(10.seconds)) should be(4)
    }

    it should "clean up & close the database" taggedAs DbmsTest in {
      // Not relevant in Travis where all state gets nuked but useful for testing locally
      database.runTestTransaction(database.dataAccess.metadataEntries.delete).futureValue(Timeout(10.seconds))
      database.runTestTransaction(database.dataAccess.workflowMetadataSummaryEntries.delete).futureValue(Timeout(10.seconds))

      database.close()
    }

    it should "stop container if required" taggedAs DbmsTest in {
      containerOpt.foreach { _.stop }
    }

  }

  behavior of "MetadataSlickDatabase"
  it should "partition metadata for summarization correctly" in {

    def partition(metadata: Seq[MetadataEntry]): SummarizationPartitionedMetadata = {
      MetadataSlickDatabase.partitionSummarizationMetadata(
        rawMetadataEntries = metadata,
        startMetadataKey = WorkflowMetadataKeys.StartTime,
        endMetadataKey = WorkflowMetadataKeys.EndTime,
        nameMetadataKey = WorkflowMetadataKeys.Name,
        statusMetadataKey = WorkflowMetadataKeys.Status,
        submissionMetadataKey = WorkflowMetadataKeys.SubmissionTime,
        parentWorkflowIdKey = WorkflowMetadataKeys.ParentWorkflowId,
        rootWorkflowIdKey = WorkflowMetadataKeys.RootWorkflowId,
        labelMetadataKey = WorkflowMetadataKeys.Labels)
    }

    {
      // Edge condition: empty input
      val partitioned = partition(List.empty)
      partitioned.nonSummarizableMetadata shouldBe empty
      partitioned.summarizableMetadata shouldBe empty
    }

    {
      // A mix of summarizable and non-summarizable keys specified at workflow and call levels.
      val wfid = WorkflowId.randomId().id.toString
      val callName = "my.call"

      def callEntry(key: String): MetadataEntry =
        MetadataEntry(wfid, Option(callName), None, Option(1), key, None, None, OffsetDateTime.now().toSystemTimestamp)

      def workflowEntry(key: String): MetadataEntry =
        MetadataEntry(wfid, None, None, None, key, None, None, OffsetDateTime.now().toSystemTimestamp)

      val rightKeysCallLevel = List(
        callEntry(WorkflowMetadataKeys.StartTime),
        callEntry(WorkflowMetadataKeys.EndTime),
        callEntry(WorkflowMetadataKeys.Name),
        callEntry(WorkflowMetadataKeys.Status),
        callEntry(WorkflowMetadataKeys.SubmissionTime),
        callEntry(WorkflowMetadataKeys.ParentWorkflowId),
        callEntry(WorkflowMetadataKeys.RootWorkflowId),
        callEntry(WorkflowMetadataKeys.Labels + ":arbitrary-label")
      )

      val wrongKeysCallLevel = List(
        callEntry("complete"),
        callEntry("rubbish")
      )

      val thingsThatLookKindOfLikeTheRightWorkflowKeysButActuallyAreNotAndAreCallScopedAnyway = List(
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.StartTime),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.EndTime),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.Name),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.Status),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.SubmissionTime),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.ParentWorkflowId),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.RootWorkflowId),
        callEntry(CallMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.Labels + ":arbitrary-label"),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.StartTime),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.EndTime),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.Name),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.Status),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.SubmissionTime),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.ParentWorkflowId),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.RootWorkflowId),
        callEntry(CallMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.Labels + ":arbitrary-label")
      )

      val rightKeysWorkflowLevel = List(
        workflowEntry(WorkflowMetadataKeys.StartTime),
        workflowEntry(WorkflowMetadataKeys.EndTime),
        workflowEntry(WorkflowMetadataKeys.Name),
        workflowEntry(WorkflowMetadataKeys.Status),
        workflowEntry(WorkflowMetadataKeys.SubmissionTime),
        workflowEntry(WorkflowMetadataKeys.ParentWorkflowId),
        workflowEntry(WorkflowMetadataKeys.RootWorkflowId),
        workflowEntry(WorkflowMetadataKeys.Labels + ":arbitrary-label")
      )

      val wrongKeysWorkflowLevel = List(
        workflowEntry("total"),
        workflowEntry("garbage")
      )

      val thingsThatLookKindOfLikeTheRightWorkflowKeysButActuallyAreNot = List(
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.StartTime),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.EndTime),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.Name),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.Status),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.SubmissionTime),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.ParentWorkflowId),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.RootWorkflowId),
        workflowEntry(WorkflowMetadataKeys.Inputs + ":" + WorkflowMetadataKeys.Labels + ":arbitrary-label"),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.StartTime),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.EndTime),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.Name),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.Status),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.SubmissionTime),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.ParentWorkflowId),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.RootWorkflowId),
        workflowEntry(WorkflowMetadataKeys.Outputs + ":" + WorkflowMetadataKeys.Labels + ":arbitrary-label")
      )

      val allTheWrongThings = rightKeysCallLevel ++ wrongKeysCallLevel ++ wrongKeysWorkflowLevel ++
        thingsThatLookKindOfLikeTheRightWorkflowKeysButActuallyAreNot ++
        thingsThatLookKindOfLikeTheRightWorkflowKeysButActuallyAreNotAndAreCallScopedAnyway

      val partitioned = partition(rightKeysWorkflowLevel ++ allTheWrongThings)
      partitioned.nonSummarizableMetadata.toSet shouldBe (allTheWrongThings).toSet
      partitioned.summarizableMetadata shouldBe rightKeysWorkflowLevel
    }
  }
}
