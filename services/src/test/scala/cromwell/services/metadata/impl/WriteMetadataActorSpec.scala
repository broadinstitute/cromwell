package cromwell.services.metadata.impl

import java.sql.{Connection, Timestamp}

import akka.testkit.{TestFSMRef, TestProbe}
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.database.sql.joins.MetadataJobQueryValue
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.database.sql.{MetadataSqlDatabase, SqlDatabase, StreamMetadataSqlDatabase}
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
  val mockDatabaseInterface = new MetadataSqlDatabase with SqlDatabase with StreamMetadataSqlDatabase {
    private def notImplemented() = throw new UnsupportedOperationException

    override protected val urlKey = "mock_database_url"
    override protected val originalDatabaseConfig = ConfigFactory.empty

    override def connectionDescription: String = "Mock Database"

    override def existsMetadataEntries()(
      implicit ec: ExecutionContext): Nothing = notImplemented()

    // Return successful
    override def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])
                                   (implicit ec: ExecutionContext): Future[Unit] = Future.successful(())

    override def metadataEntryExists(workflowExecutionUuid: String)
                                    (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntries(workflowExecutionUuid: String)
                                     (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntries(workflowExecutionUuid: String,
                                      metadataKey: String)(implicit ec: ExecutionContext): Nothing = {
      notImplemented()
    }

    override def queryMetadataEntries(workflowExecutionUuid: String,
                                      callFullyQualifiedName: String,
                                      jobIndex: Option[Int],
                                      jobAttempt: Option[Int])(implicit ec: ExecutionContext): Nothing = {
      notImplemented()
    }

    override def queryMetadataEntries(workflowUuid: String,
                                      metadataKey: String,
                                      callFullyQualifiedName: String,
                                      jobIndex: Option[Int],
                                      jobAttempt: Option[Int])(implicit ec: ExecutionContext): Nothing = {
      notImplemented()
    }

    override def queryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                                      metadataKeys: NonEmptyList[String],
                                                      metadataJobQueryValue: MetadataJobQueryValue)
                                                     (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                                       metadataKeys: NonEmptyList[String],
                                                       metadataJobQueryValue: MetadataJobQueryValue)
                                                      (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def refreshMetadataSummaryEntries(startMetadataKey: String,
                                               endMetadataKey: String,
                                               nameMetadataKey: String,
                                               statusMetadataKey: String,
                                               labelMetadataKey: String,
                                               submissionMetadataKey: String,
                                               buildUpdatedSummary: (
                                                 Option[WorkflowMetadataSummaryEntry],
                                                   Seq[MetadataEntry]) => WorkflowMetadataSummaryEntry)
                                              (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def getWorkflowStatus(workflowExecutionUuid: String)
                                  (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def getWorkflowLabels(workflowExecutionUuid: String)
                                  (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryWorkflowSummaries(parentWorkflowIdMetadataKey: String,
                                        workflowStatuses: Set[String],
                                        workflowNames: Set[String],
                                        workflowExecutionUuids: Set[String],
                                        labelAndKeyLabelValues: Set[(String, String)],
                                        labelOrKeyLabelValues: Set[(String, String)],
                                        excludeLabelAndValues: Set[(String, String)],
                                        excludeLabelOrValues: Set[(String, String)],
                                        submissionTimestamp: Option[Timestamp],
                                        startTimestampOption: Option[Timestamp],
                                        endTimestampOption: Option[Timestamp],
                                        includeSubworkflows: Boolean,
                                        page: Option[Int],
                                        pageSize: Option[Int])
                                       (implicit ec: ExecutionContext): Nothing = {
      notImplemented()
    }

    override def countWorkflowSummaries(parentWorkflowIdMetadataKey: String,
                                        workflowStatuses: Set[String],
                                        workflowNames: Set[String],
                                        workflowExecutionUuids: Set[String],
                                        labelAndKeyLabelValues: Set[(String, String)],
                                        labelOrKeyLabelValues: Set[(String, String)],
                                        excludeLabelAndValues: Set[(String, String)],
                                        excludeLabelOrValues: Set[(String, String)],
                                        submissionTimestamp: Option[Timestamp],
                                        startTimestampOption: Option[Timestamp],
                                        endTimestampOption: Option[Timestamp],
                                        includeSubworkflows: Boolean)(implicit ec: ExecutionContext): Nothing = {
      notImplemented()
    }

    override def withConnection[A](block: Connection => A): Nothing = {
      notImplemented()
    }

    override def close(): Nothing = notImplemented()

    override def streamQueryMetadataEntries(workflowExecutionUuid: String): Nothing = notImplemented()

    override def streamQueryMetadataEntries(workflowExecutionUuid: String, metadataKey: String): Nothing = notImplemented()

    override def streamQueryMetadataEntries(workflowExecutionUuid: String, callFullyQualifiedName: String, jobIndex: Option[Int], jobAttempt: Option[Int]): Nothing = notImplemented()

    override def streamQueryMetadataEntries(workflowUuid: String, metadataKey: String, callFullyQualifiedName: String, jobIndex: Option[Int], jobAttempt: Option[Int]): Nothing = notImplemented()

    override def streamQueryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String], metadataJobQueryValue: MetadataJobQueryValue): Nothing = notImplemented()

    override def streamQueryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String, metadataKeys: NonEmptyList[String], metadataJobQueryValue: MetadataJobQueryValue): Nothing = notImplemented()
  }
}
