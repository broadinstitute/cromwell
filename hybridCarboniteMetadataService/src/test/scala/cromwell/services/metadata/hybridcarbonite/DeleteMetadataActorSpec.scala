package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime
import java.util.UUID

import akka.actor.Props
import akka.testkit.{EventFilter, ImplicitSender, TestProbe}
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus.Archived
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor.DeleteMetadataAction
import org.scalatest.flatspec.AnyFlatSpecLike

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

class DeleteMetadataActorSpec extends TestKitSuite with AnyFlatSpecLike with ImplicitSender {

  private val workflowId1 = UUID.randomUUID().toString
  private val workflowId2 = UUID.randomUUID().toString

  val metricActor = TestProbe().ref

  it should "lookup for workflows matching criteria and delete their metadata" in {
    val deleteMetadataActor = createTestDeletionActor()
    EventFilter.info(pattern = s"Successfully deleted metadata for workflow *", occurrences = 2) intercept {
      deleteMetadataActor ! DeleteMetadataAction
    }
  }

  it should "write error message to the log if lookup fails" in {
    val deleteMetadataActorFailingLookups = createTestDeletionActor(failLookups = true)
    EventFilter.error(start = "Cannot delete metadata: unable to query list of workflow ids for metadata deletion from metadata summary table.", occurrences = 1) intercept {
      deleteMetadataActorFailingLookups ! DeleteMetadataAction
    }
  }

  it should "write error message to the log if deletion fails and continue processing next workflow from the list" in {
    val deleteMetadataActorFailingDeletions = createTestDeletionActor(failDeletions = true)
    // error message should appear twice
    EventFilter.error(pattern = s"Cannot delete metadata for workflow *", occurrences = 2) intercept {
      deleteMetadataActorFailingDeletions ! DeleteMetadataAction
    }
  }

  private def createTestDeletionActor(failLookups: Boolean = false, failDeletions: Boolean = false) = {
    val deleteMetadataActor = system.actorOf(Props(new DeleteMetadataActor(ActiveMetadataDeletionConfig(1 minute, 1 minute, 200L, 1 minute), serviceRegistryActor = TestProbe().ref) {
      override def queryRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(archiveStatus: Option[String], thresholdTimestamp: OffsetDateTime, batchSize: Long)(implicit ec: ExecutionContext): Future[Seq[String]] = {
        val expectedArchiveStatus = MetadataArchiveStatus.toDatabaseValue(Archived)
        if (expectedArchiveStatus !== archiveStatus) {
          Future.failed(new RuntimeException(s"DeleteMetadataActor misconfiguration: should query entries with $expectedArchiveStatus status, instead queries $archiveStatus"))
        } else if (thresholdTimestamp.isAfter(OffsetDateTime.now().minusMinutes(1))) {
          Future.failed(new RuntimeException(s"DeleteMetadataActor misconfiguration: thresholdTimestamp should be set to at least a minute in the past"))
        } else if (failLookups) {
          Future.failed(new RuntimeException("Error occurred during lookup for metadata deletion candidates"))
        } else {
          Future.successful(Seq(workflowId1, workflowId2))
        }
      }

      override def deleteNonLabelMetadataEntriesForWorkflowAndUpdateArchiveStatus(rootWorkflowId: WorkflowId, newArchiveStatus: Option[String])(implicit ec: ExecutionContext): Future[Int] = {
        if (failDeletions) {
          Future.failed(new RuntimeException(s"Error occurred during metadata deletion for workflow $rootWorkflowId"))
        } else {
          Future.successful(1)
        }
      }
    }))
    deleteMetadataActor
  }
}


