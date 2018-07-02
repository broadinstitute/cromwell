package cromwell.services.metadata.impl

import java.sql.{Connection, Timestamp}

import akka.testkit.{TestFSMRef, TestProbe}
import cats.data.NonEmptyList
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.database.sql.joins.MetadataJobQueryValue
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.database.sql.{MetadataSqlDatabase, SqlDatabase}
import cromwell.services.metadata.MetadataService.{MetadataWriteSuccess, PutMetadataActionAndRespond}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

class WriteMetadataActorSpec extends TestKitSuite with FlatSpecLike with Matchers {

  behavior of "WriteMetadataActor"
  
  it should "respond to all senders" in {
    val registry = TestProbe().ref
    val writeActor = TestFSMRef(new WriteMetadataActor(10, 1.second, registry, Int.MaxValue) {
      override val metadataDatabaseInterface = mockDatabaseInterface
    })
    
    val metadataKey = MetadataKey(WorkflowId.randomId(), None, "metadata_key") 
    val metadataEvent = MetadataEvent(metadataKey, MetadataValue("hello"))
    
    val probes = (1 until 20).map({ _ =>
      val probe = TestProbe()
      probe.send(writeActor, PutMetadataActionAndRespond(List(metadataEvent), probe.ref))
      probe
    })
    
    probes.foreach(_.expectMsg(MetadataWriteSuccess(List(metadataEvent))))
  }

  // Mock database interface.
  val mockDatabaseInterface = new MetadataSqlDatabase with SqlDatabase {
    // Return successful
    override def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])(implicit ec: ExecutionContext) = Future.successful(())

    override def metadataEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext) = ???
    override def queryMetadataEntries(workflowExecutionUuid: String)(implicit ec: ExecutionContext) = ???
    override def queryMetadataEntries(workflowExecutionUuid: String, metadataKey: String)(implicit ec: ExecutionContext) = ???
    override def queryMetadataEntries(workflowExecutionUuid: String, callFullyQualifiedName: String, jobIndex: Option[Int], jobAttempt: Option[Int])(implicit ec: ExecutionContext) = ???
    override def queryMetadataEntries(workflowUuid: String, metadataKey: String, callFullyQualifiedName: String, jobIndex: Option[Int], jobAttempt: Option[Int])(implicit ec: ExecutionContext) = ???
    override def queryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String], metadataJobQueryValue: MetadataJobQueryValue)(implicit ec: ExecutionContext) = ???
    override def queryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String], metadataJobQueryValue: MetadataJobQueryValue)(implicit ec: ExecutionContext) = ???
    override def refreshMetadataSummaryEntries(startMetadataKey: String, endMetadataKey: String, nameMetadataKey: String, statusMetadataKey: String, labelMetadataKey: String, submissionMetadataKey: String, buildUpdatedSummary: (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry]) => WorkflowMetadataSummaryEntry)(implicit ec: ExecutionContext) = ???
    override def getWorkflowStatus(workflowExecutionUuid: String)(implicit ec: ExecutionContext) = ???
    override def getWorkflowLabels(workflowExecutionUuid: String)(implicit ec: ExecutionContext) = ???
    override def queryWorkflowSummaries(workflowStatuses: Set[String], workflowNames: Set[String], workflowExecutionUuids: Set[String], labelAndKeyLabelValues: Set[(String, String)], labelOrKeyLabelValues: Set[(String, String)], excludeLabelAndValues: Set[(String, String)], excludeLabelOrValues: Set[(String, String)], submissionTimestamp: Option[Timestamp], startTimestampOption: Option[Timestamp], endTimestampOption: Option[Timestamp], page: Option[Int], pageSize: Option[Int])(implicit ec: ExecutionContext) = ???
    override def countWorkflowSummaries(workflowStatuses: Set[String], workflowNames: Set[String], workflowExecutionUuids: Set[String], labelAndKeyLabelValues: Set[(String, String)], labelOrKeyLabelValues: Set[(String, String)], excludeLabelAndValues: Set[(String, String)], excludeLabelOrValues: Set[(String, String)], submissionTimestamp: Option[Timestamp], startTimestampOption: Option[Timestamp], endTimestampOption: Option[Timestamp])(implicit ec: ExecutionContext) = ???
    override protected val urlKey = null
    override protected val originalDatabaseConfig = null
    override def withConnection[A](block: Connection => A) = ???
    override def close() = ???
  }

}
