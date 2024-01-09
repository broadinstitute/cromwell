package cromwell.services.database

import java.time.OffsetDateTime
import com.dimafeng.testcontainers.Container
import common.assertion.CromwellTimeoutSpec
import cromwell.core.Tags.DbmsTest
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.database.migration.metadata.table.symbol.MetadataStatement._
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.database.slick.MetadataSlickDatabase.SummarizationPartitionedMetadata
import cromwell.database.sql.joins.{CallOrWorkflowQuery, CallQuery, WorkflowQuery}
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.services.metadata.CallMetadataKeys
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import javax.sql.rowset.serial.SerialClob
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future}
import scala.language.postfixOps

class MetadataSlickDatabaseSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with ScalaFutures {

  DatabaseSystem.All foreach { databaseSystem =>
    implicit val ec: ExecutionContextExecutor = ExecutionContext.global

    behavior of s"MetadataSlickDatabase on ${databaseSystem.name}"

    val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

    lazy val database =
      DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, MetadataDatabaseType, databaseSystem)
    import cromwell.database.migration.metadata.table.symbol.MetadataStatement.OffsetDateTimeToSystemTimestamp
    import database.dataAccess.driver.api._
    val now = OffsetDateTime.now().toSystemTimestamp

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.start)
    }

    val rootCountableId = "root workflow id: countable stuff"

    val subWorkflowCountableId = "subworkflow id: countable stuff"

    val subSubWorkflowCountableId = "subsubworkflow id: countable stuff"

    val failedParentWorkflowId = "bbf4c25b-282b-4a18-a914-441f9684b69e"
    val ignoredFailedParentWorkflowId = "ccccc-ddddd-3333-42242"
    val ignoredFailedChildWorkflowId = "eeeeeee-ggggggg-33333-22222"
    val successfulParentWorkflowId = "4c1cf43d-1fbd-47af-944c-c63216f293ae"

    val failedChildWorkflowId = "9ff3d855-0585-48e4-b3a1-189101f611e5"
    val successfulChildWorkflowId = "73886096-2e06-48f6-ba42-f365dbf23de5"

    val failedStatusMetadataValue = Option(new SerialClob("Failed".toCharArray()))
    val doneStatusMetadataValue = Option(new SerialClob("Done".toCharArray()))
    val stdErrValue = Option(new SerialClob("test value".toCharArray()))

    it should "set up the test data" taggedAs DbmsTest in {
      database
        .runTestTransaction(
          database.dataAccess.metadataEntries ++= Seq(
            MetadataEntry("workflow id: 4 to delete, including 1 label",
                          None,
                          None,
                          None,
                          "someKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: 4 to delete, including 1 label",
                          None,
                          None,
                          None,
                          "someKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: 4 to delete, including 1 label",
                          None,
                          None,
                          None,
                          "someKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: 4 to delete, including 1 label",
                          None,
                          None,
                          None,
                          "labels:do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: I am a root workflow with a subworkflow",
                          None,
                          None,
                          None,
                          "labels:do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: I am a root workflow with a subworkflow",
                          None,
                          None,
                          None,
                          "please do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: I am the subworkflow",
                          None,
                          None,
                          None,
                          "labels:do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("workflow id: I am the subworkflow",
                          None,
                          None,
                          None,
                          "please do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("nested subworkflows: root",
                          None,
                          None,
                          None,
                          "please do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("nested subworkflows: first nesting",
                          None,
                          None,
                          None,
                          "please do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("nested subworkflows: second nesting",
                          None,
                          None,
                          None,
                          "please do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry("nested subworkflows: third nesting",
                          None,
                          None,
                          None,
                          "please do delete me",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),

            // Workflow level
            MetadataEntry(rootCountableId,
                          None,
                          None,
                          None,
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(rootCountableId,
                          None,
                          None,
                          None,
                          "excludableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            // Call level
            MetadataEntry(rootCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(1),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(rootCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(1),
                          "excludableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            //   Excludable call index
            MetadataEntry(rootCountableId,
                          Option("includableCall"),
                          Option(1),
                          Option(1),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            //   Excludable call attempt
            MetadataEntry(rootCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(2),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(rootCountableId,
                          Option("excludableCall"),
                          Option(0),
                          Option(1),
                          "whateverKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),

            // Subworkflow level
            MetadataEntry(subWorkflowCountableId,
                          None,
                          None,
                          None,
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(subWorkflowCountableId,
                          None,
                          None,
                          None,
                          "excludableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            // Call level
            MetadataEntry(subWorkflowCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(1),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(subWorkflowCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(1),
                          "excludableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            //   Excludable call index
            MetadataEntry(subWorkflowCountableId,
                          Option("includableCall"),
                          Option(1),
                          Option(1),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            //   Excludable call attempt
            MetadataEntry(subWorkflowCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(2),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(subWorkflowCountableId,
                          Option("excludableCall"),
                          Option(0),
                          Option(1),
                          "whateverKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),

            // Subsubworkflow level
            MetadataEntry(subSubWorkflowCountableId,
                          None,
                          None,
                          None,
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(subSubWorkflowCountableId,
                          None,
                          None,
                          None,
                          "excludableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            // Call level
            MetadataEntry(subSubWorkflowCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(1),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(subSubWorkflowCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(1),
                          "excludableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            //   Excludable call index
            MetadataEntry(subSubWorkflowCountableId,
                          Option("includableCall"),
                          Option(1),
                          Option(1),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            //   Excludable call attempt
            MetadataEntry(subSubWorkflowCountableId,
                          Option("includableCall"),
                          Option(0),
                          Option(2),
                          "includableKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            ),
            MetadataEntry(subSubWorkflowCountableId,
                          Option("excludableCall"),
                          Option(0),
                          Option(1),
                          "whateverKey",
                          None,
                          None,
                          OffsetDateTime.now().toSystemTimestamp,
                          None
            )
          )
        )
        .futureValue(Timeout(10.seconds))

      database
        .runTestTransaction(
          database.dataAccess.workflowMetadataSummaryEntries ++= Seq(
            WorkflowMetadataSummaryEntry("workflow id: 3 to delete, 1 label",
                                         Option("workflow name"),
                                         Option("Succeeded"),
                                         Option(now),
                                         Option(now),
                                         Option(now),
                                         None,
                                         None,
                                         None,
                                         None
            ),
            WorkflowMetadataSummaryEntry(
              "workflow id: I am a root workflow with a subworkflow",
              Option("workflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              None,
              None,
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(
              "workflow id: I am the subworkflow",
              Option("workflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option("workflow id: I am a root workflow with a subworkflow"),
              Option("workflow id: I am a root workflow with a subworkflow"),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry("nested subworkflows: root",
                                         Option("workflow name"),
                                         Option("Succeeded"),
                                         Option(now),
                                         Option(now),
                                         Option(now),
                                         None,
                                         None,
                                         None,
                                         None
            ),
            WorkflowMetadataSummaryEntry(
              "nested subworkflows: first nesting",
              Option("workflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option("nested subworkflows: root"),
              Option("nested subworkflows: root"),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(
              "nested subworkflows: second nesting",
              Option("workflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option("nested subworkflows: first nesting"),
              Option("nested subworkflows: root"),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(
              "nested subworkflows: third nesting",
              Option("workflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option("nested subworkflows: second nesting"),
              Option("nested subworkflows: root"),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(rootCountableId,
                                         Option("workflow name"),
                                         Option("Succeeded"),
                                         Option(now),
                                         Option(now),
                                         Option(now),
                                         None,
                                         None,
                                         None,
                                         None
            ),
            WorkflowMetadataSummaryEntry(
              subWorkflowCountableId,
              Option("subworkflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option(rootCountableId),
              Option(rootCountableId),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(
              subSubWorkflowCountableId,
              Option("subsubworkflow name"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option(subWorkflowCountableId),
              Option(rootCountableId),
              None,
              None
            )
          )
        )
        .futureValue(Timeout(10.seconds))
    }

    it should "delete the right number of rows for a root workflow without subworkflows" taggedAs DbmsTest in {
      val delete =
        database.deleteAllMetadataForWorkflowAndUpdateArchiveStatus("workflow id: 4 to delete, including 1 label", None)
      delete.futureValue(Timeout(10.seconds)) should be(4)
    }

    it should "delete the right number of rows for a root workflow with subworkflows" taggedAs DbmsTest in {
      val delete = database.deleteAllMetadataForWorkflowAndUpdateArchiveStatus(
        "workflow id: I am a root workflow with a subworkflow",
        None
      )
      delete.futureValue(Timeout(10.seconds)) should be(2)
    }

    it should "delete the right number of rows for a nested subworkflow" taggedAs DbmsTest in {
      val delete = database.deleteAllMetadataForWorkflowAndUpdateArchiveStatus("nested subworkflows: root", None)
      delete.futureValue(Timeout(10.seconds)) should be(1)
    }

    it should "count up rows" taggedAs DbmsTest in {
      List(true, false) foreach { expandSubWorkflows =>
        val expansionFactor = if (expandSubWorkflows) 3 else 1;
        // Everything
        {
          val count =
            database.countMetadataEntries(rootCountableId, expandSubWorkflows = expandSubWorkflows, 10 seconds)
          count.futureValue(Timeout(10.seconds)) should be(7 * expansionFactor)
        }

        // Only includable keys - this looks for workflow level data only
        {
          val count = database.countMetadataEntries(rootCountableId,
                                                    "includableKey",
                                                    expandSubWorkflows = expandSubWorkflows,
                                                    10 seconds
          )
          count.futureValue(Timeout(10.seconds)) should be(1 * expansionFactor)
        }

        {
          val count = database.countMetadataEntries(rootCountableId,
                                                    "includableCall",
                                                    Option(0),
                                                    Option(1),
                                                    expandSubWorkflows = expandSubWorkflows,
                                                    10 seconds
          )
          count.futureValue(Timeout(10 seconds)) should be(2 * expansionFactor)
        }

        {
          val count = database.countMetadataEntries(rootCountableId,
                                                    "includableKey",
                                                    "includableCall",
                                                    Option(0),
                                                    Option(1),
                                                    expandSubWorkflows = expandSubWorkflows,
                                                    10 seconds
          )
          count.futureValue(Timeout(10 seconds)) should be(1 * expansionFactor)
        }

        {
          val count = database.countMetadataEntryWithKeyConstraints(
            workflowExecutionUuid = rootCountableId,
            metadataKeysToFilterFor = List("includable%"),
            metadataKeysToFilterOut = List("excludable%"),
            CallQuery("includableCall", Option(0), Option(1)),
            expandSubWorkflows = expandSubWorkflows,
            10 seconds
          )
          count.futureValue(Timeout(10 seconds)) should be(1 * expansionFactor)
        }

        {
          val count = database.countMetadataEntryWithKeyConstraints(
            workflowExecutionUuid = rootCountableId,
            metadataKeysToFilterFor = List("includable%"),
            metadataKeysToFilterOut = List("excludable%"),
            CallOrWorkflowQuery,
            expandSubWorkflows = expandSubWorkflows,
            10 seconds
          )
          count.futureValue(Timeout(10 seconds)) should be(4 * expansionFactor)
        }

        {
          val count = database.countMetadataEntryWithKeyConstraints(
            workflowExecutionUuid = rootCountableId,
            metadataKeysToFilterFor = List("includable%"),
            metadataKeysToFilterOut = List("excludable%"),
            WorkflowQuery,
            expandSubWorkflows = expandSubWorkflows,
            10 seconds
          )
          count.futureValue(Timeout(10 seconds)) should be(1 * expansionFactor)
        }
      }
    }

    it should "fetch failed tasks from a failed workflow" taggedAs DbmsTest in {
      database
        .runTestTransaction(
          database.dataAccess.metadataEntries ++= Seq(
            // Failed parent workflow calls (successful calls and older attempts and runs mixed in for negative checks)
            MetadataEntry(
              failedParentWorkflowId,
              Option("failedWorkflowCall"),
              Option(0),
              Option(0),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedParentWorkflowId,
              Option("failedWorkflowCall"),
              Option(0),
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedParentWorkflowId,
              Option("failedWorkflowCall"),
              Option(1),
              Option(0),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedParentWorkflowId,
              Option("failedWorkflowCall"),
              Option(1),
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedParentWorkflowId,
              Option("failedWorkflowCall"),
              Option(1),
              Option(1),
              "stderr",
              stdErrValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedParentWorkflowId,
              Option("successfulWorkflowCall"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedParentWorkflowId,
              Option("successfulSubWorkflowCall"),
              Option(0),
              Option(0),
              "subWorkflowId",
              Option(new SerialClob(failedChildWorkflowId.toCharArray)),
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // ignored failed workflow calls. These should not be fetched since it's not part of the target workflow tree
            MetadataEntry(
              ignoredFailedParentWorkflowId,
              Option("failedWorkflowCall"),
              Option(0),
              Option(0),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              ignoredFailedChildWorkflowId,
              Option("failedSubWorkflowCall"),
              None,
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // child workflow calls (successful calls and previous failed attempts/shards are mixed in for negative checks)
            // notFailedSubWorkflowCall should not be returned since it succeeded on the last attempt and has no scatters
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallyFailedSubWorkflowCall"),
              None,
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallyFailedSubWorkflowCall"),
              None,
              Option(2),
              "backendStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallyFailedSubWorkflowCall"),
              None,
              Option(2),
              "stderr",
              stdErrValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("successfulSubWorkflowCall"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // Failed child workflow calls (successful calls and previous failed attempts/shards are mixed in for negative checks)
            // failedSubWorkflowCall should be returned since it never succeeded
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall"),
              None,
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall"),
              None,
              Option(2),
              "backendStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall"),
              None,
              Option(2),
              "stderr",
              stdErrValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("successfulSubWorkflowCall2"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // Third set of child workflow calls, similar to above however this set consists of retries and scatters
            // It's safe to assume that if one scatter fails then they all fail, so pull the last scatter on the last attempt
            // failedSubWorkflowCall2 should return since the scatters failed on the last attempt
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall2"),
              Option(1),
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall2"),
              Option(2),
              Option(1),
              "backendStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall2"),
              Option(1),
              Option(2),
              "backendStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall2"),
              Option(2),
              Option(2),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("failedSubWorkflowCall2"),
              Option(2),
              Option(2),
              "stderr",
              stdErrValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("successfulSubWorkflowCall3"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // Fourth set of child workflow calls
            // This set should not return anything since the scatters succeeded on the last attempt
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallySubWorkflowCall2"),
              Option(1),
              Option(1),
              "executionStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallySubWorkflowCall2"),
              Option(2),
              Option(1),
              "backendStatus",
              failedStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallySubWorkflowCall2"),
              Option(1),
              Option(2),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallySubWorkflowCall2"),
              Option(2),
              Option(2),
              "backendStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("notActuallySubWorkflowCall2"),
              Option(2),
              Option(2),
              "stderr",
              stdErrValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              failedChildWorkflowId,
              Option("successfulSubWorkflowCall4"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // Successful parent workflow call (negative check)
            MetadataEntry(
              successfulParentWorkflowId,
              Option("successfulWorkflowCall"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),

            // Successful child workflow calls (negative check)
            MetadataEntry(
              successfulChildWorkflowId,
              Option("successfulSubWorkflowCall"),
              Option(0),
              Option(0),
              "executionStatus",
              doneStatusMetadataValue,
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            ),
            MetadataEntry(
              successfulChildWorkflowId,
              Option("successfulSubWorkflowCall"),
              Option(0),
              Option(0),
              "subWorkflowId",
              Option(new SerialClob(successfulParentWorkflowId.toCharArray)),
              Option("String"),
              OffsetDateTime.now().toSystemTimestamp,
              None
            )
          )
        )
        .futureValue(Timeout(10.seconds))

      database
        .runTestTransaction(
          database.dataAccess.workflowMetadataSummaryEntries ++= Seq(
            // Failed WorkflowMetadataSummaryEntry setup
            WorkflowMetadataSummaryEntry(failedParentWorkflowId,
                                         Option("failedParentWorkflow"),
                                         Option("Failed"),
                                         Option(now),
                                         Option(now),
                                         Option(now),
                                         None,
                                         None,
                                         None,
                                         None
            ),
            WorkflowMetadataSummaryEntry(
              failedChildWorkflowId,
              Option("failedChildWorkflow"),
              Option("Failed"),
              Option(now),
              Option(now),
              Option(now),
              Option(failedParentWorkflowId),
              Option(failedParentWorkflowId),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(successfulParentWorkflowId,
                                         Option("successfulParentWorkflow"),
                                         Option("Succeeded"),
                                         Option(now),
                                         Option(now),
                                         Option(now),
                                         None,
                                         None,
                                         None,
                                         None
            ),
            WorkflowMetadataSummaryEntry(
              successfulChildWorkflowId,
              Option("successfulChildWorkflow"),
              Option("Succeeded"),
              Option(now),
              Option(now),
              Option(now),
              Option(successfulParentWorkflowId),
              Option(successfulParentWorkflowId),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(
              ignoredFailedParentWorkflowId,
              Option("ignoredFailedParentWorkflow"),
              Option("Failed"),
              Option(now),
              Option(now),
              Option(now),
              Option(ignoredFailedParentWorkflowId),
              Option(ignoredFailedParentWorkflowId),
              None,
              None
            ),
            WorkflowMetadataSummaryEntry(
              ignoredFailedChildWorkflowId,
              Option("ignoredFailedChildWorkflow"),
              Option("Failed"),
              Option(now),
              Option(now),
              Option(now),
              Option(ignoredFailedChildWorkflowId),
              Option(ignoredFailedChildWorkflowId),
              None,
              None
            )
          )
        )
        .futureValue(Timeout(10.seconds))

      val futureEntries: Future[Seq[MetadataEntry]] =
        database.getFailedJobsMetadataWithWorkflowId(failedParentWorkflowId)
      var entriesFound = false

      val recordCount = Map(
        failedParentWorkflowId -> scala.collection.mutable.Map(
          "expected" -> 2,
          "actual" -> 0
        ),
        failedChildWorkflowId -> scala.collection.mutable.Map(
          "expected" -> 4,
          "actual" -> 0
        )
      )

      whenReady(futureEntries, timeout(scaled(5 seconds))) { entries =>
        entries.foreach { entry =>
          entriesFound = true
          val workflowId = entry.workflowExecutionUuid
          recordCount.getOrElse(workflowId, None) should not be None
          recordCount(workflowId)("actual") += 1
          val metadataValueClob = entry.metadataValue.get
          val metadataValueString = metadataValueClob.getSubString(1, metadataValueClob.length().toInt)

          entry.metadataKey match {
            case "stderr" =>
              metadataValueString should be("test value")
            case "backendStatus" =>
              metadataValueString should be("Failed")
            case "executionStatus" =>
              metadataValueString should be("Failed")
            case _ => fail(s"Unexpected key ${entry.metadataKey} was included in result set")
          }

          entry.metadataKey should not be "subWorkflowId"
          entry.callFullyQualifiedName.getOrElse("") match {
            case "failedWorkflowCall" =>
              entry.jobIndex.get should be(1)
              entry.jobAttempt.get should be(1)
              val isValidKey = List("executionStatus", "stderr").contains(entry.metadataKey)
              isValidKey should be(true)
            case "failedSubWorkflowCall" =>
              entry.jobIndex should be(None)
              entry.jobAttempt.get should be(2)
              val isValidKey = List("backendStatus", "stderr").contains(entry.metadataKey)
              isValidKey should be(true)
            case "failedSubWorkflowCall2" =>
              entry.jobIndex.get should be(2)
              entry.jobAttempt.get should be(2)
              val isValidKey = List("executionStatus", "stderr").contains(entry.metadataKey)
              isValidKey should be(true)
            case _ =>
              fail(
                s"Entry ${entry.callFullyQualifiedName.getOrElse("N/A")} | Index: ${entry.jobIndex.get} | Attempt: ${entry.jobAttempt.get} should not be in result set"
              )
          }
        }
        entriesFound should be(true)
        recordCount.foreach(record => record._2("actual") should be(record._2("expected")))
      }
    }

    it should "clean up & close the database" taggedAs DbmsTest in {
      // Not relevant in Travis where all state gets nuked but useful for testing locally
      database.runTestTransaction(database.dataAccess.metadataEntries.delete).futureValue(Timeout(10.seconds))
      database
        .runTestTransaction(database.dataAccess.workflowMetadataSummaryEntries.delete)
        .futureValue(Timeout(10.seconds))

      database.close()
    }

    it should "stop container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.stop())
    }

  }

  behavior of "MetadataSlickDatabase"
  it should "partition metadata for summarization correctly" in {

    def partition(metadata: Seq[MetadataEntry]): SummarizationPartitionedMetadata =
      MetadataSlickDatabase.partitionSummarizationMetadata(
        rawMetadataEntries = metadata,
        startMetadataKey = WorkflowMetadataKeys.StartTime,
        endMetadataKey = WorkflowMetadataKeys.EndTime,
        nameMetadataKey = WorkflowMetadataKeys.Name,
        statusMetadataKey = WorkflowMetadataKeys.Status,
        submissionMetadataKey = WorkflowMetadataKeys.SubmissionTime,
        parentWorkflowIdKey = WorkflowMetadataKeys.ParentWorkflowId,
        rootWorkflowIdKey = WorkflowMetadataKeys.RootWorkflowId,
        labelMetadataKey = WorkflowMetadataKeys.Labels
      )

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
      partitioned.nonSummarizableMetadata.toSet shouldBe allTheWrongThings.toSet
      partitioned.summarizableMetadata shouldBe rightKeysWorkflowLevel
    }
  }
}
