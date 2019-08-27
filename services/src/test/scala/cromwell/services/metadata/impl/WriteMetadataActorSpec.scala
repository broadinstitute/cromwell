package cromwell.services.metadata.impl

import java.sql.{Connection, Timestamp}

import akka.actor.ActorRef
import akka.testkit.{TestFSMRef, TestProbe}
import cats.data.NonEmptyVector
import com.typesafe.config.ConfigFactory
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.database.sql.joins.MetadataJobQueryValue
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.database.sql.{MetadataSqlDatabase, SqlDatabase}
import cromwell.services.metadata.MetadataService.{MetadataWriteAction, MetadataWriteFailure, MetadataWriteSuccess, PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata.impl.WriteMetadataActorSpec.BatchSizeCountingWriteMetadataActor
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import org.scalatest.concurrent.Eventually
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.control.NoStackTrace
import scala.util.{Failure, Success}

class WriteMetadataActorSpec extends TestKitSuite with FlatSpecLike with Matchers with Eventually {

  behavior of "WriteMetadataActor"

  implicit val defaultPatience = PatienceConfig(timeout = scaled(Span(10, Seconds)), interval = Span(500, Millis))

  it should "process jobs in the correct batch sizes" in {
    val registry = TestProbe().ref
    val writeActor = TestFSMRef(new BatchSizeCountingWriteMetadataActor(10, 10.millis, registry, Int.MaxValue) {
      override val metadataDatabaseInterface = mockDatabaseInterface(0)
    })

    def metadataEvent(index: Int) = PutMetadataAction(MetadataEvent(MetadataKey(WorkflowId.randomId(), None, s"metadata_key_$index"), MetadataValue(s"hello_$index")))

    val probes = (0 until 27).map({ _ =>
      val probe = TestProbe()
      probe
    }).zipWithIndex.map {
      case (probe, index) => probe -> metadataEvent(index)
    }

    probes foreach {
      case (probe, msg) => probe.send(writeActor, msg)
    }

    eventually {
      writeActor.underlyingActor.batchSizes should be(Vector(10, 10, 7))
    }

    writeActor.stop()
  }

  val failuresBetweenSuccessValues = List(0, 5, 9)
  failuresBetweenSuccessValues foreach { failureRate =>
    it should s"succeed metadata writes and respond to all senders even with $failureRate failures between each success" in {
      val registry = TestProbe().ref
      val writeActor = TestFSMRef(new BatchSizeCountingWriteMetadataActor(10, 10.millis, registry, Int.MaxValue) {
        override val metadataDatabaseInterface = mockDatabaseInterface(failureRate)
      })

      def metadataEvent(index: Int, probe: ActorRef) = PutMetadataActionAndRespond(List(MetadataEvent(MetadataKey(WorkflowId.randomId(), None, s"metadata_key_$index"), MetadataValue(s"hello_$index"))), probe)

      val probes = (0 until 43).map({ _ =>
        val probe = TestProbe()
        probe
      }).zipWithIndex.map {
        case (probe, index) => probe -> metadataEvent(index, probe.ref)
      }

      probes foreach {
        case (probe, msg) => probe.send(writeActor, msg)
      }

      probes.foreach {
        case (probe, msg) => probe.expectMsg(MetadataWriteSuccess(msg.events))
      }
      eventually {
        writeActor.underlyingActor.failureCount should be (5 * failureRate)
      }

      writeActor.stop()
    }
  }

  it should s"fail metadata writes and respond to all senders with failures" in {
    val registry = TestProbe().ref
    val writeActor = TestFSMRef(new BatchSizeCountingWriteMetadataActor(10, 10.millis, registry, Int.MaxValue) {
      override val metadataDatabaseInterface = mockDatabaseInterface(100)
    })

    def metadataEvent(index: Int, probe: ActorRef) = PutMetadataActionAndRespond(List(MetadataEvent(MetadataKey(WorkflowId.randomId(), None, s"metadata_key_$index"), MetadataValue(s"hello_$index"))), probe)

    val probes = (0 until 43).map({ _ =>
      val probe = TestProbe()
      probe
    }).zipWithIndex.map {
      case (probe, index) => probe -> metadataEvent(index, probe.ref)
    }

    probes foreach {
      case (probe, msg) => probe.send(writeActor, msg)
    }

    probes.foreach {
      case (probe, msg) => probe.expectMsg(MetadataWriteFailure(WriteMetadataActorSpec.IntermittentException, msg.events))
    }
    eventually {
      writeActor.underlyingActor.failureCount should be (5 * 10)
    }

    writeActor.stop()
  }

  // Mock database interface.
  // A customizable number of failures occur between each success
  def mockDatabaseInterface(failuresBetweenEachSuccess: Int) = new MetadataSqlDatabase with SqlDatabase {
    private def notImplemented() = throw new UnsupportedOperationException

    override protected val urlKey = "mock_database_url"
    override protected val originalDatabaseConfig = ConfigFactory.empty

    override def connectionDescription: String = "Mock Database"

    override def existsMetadataEntries()(
      implicit ec: ExecutionContext): Nothing = notImplemented()

    var requestsSinceLastSuccess = 0
    // Return successful
    override def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])
                                   (implicit ec: ExecutionContext): Future[Unit] = {
      if (requestsSinceLastSuccess == failuresBetweenEachSuccess) {
        requestsSinceLastSuccess = 0
        Future.successful(())
      } else {
        requestsSinceLastSuccess += 1
        Future.failed(WriteMetadataActorSpec.IntermittentException)
      }
    }

    override def metadataEntryExists(workflowExecutionUuid: String)
                                    (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def metadataSummaryEntryExists(workflowExecutionUuid: String)
                                           (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntries(workflowExecutionUuid: String,
                                      timeout: Duration)
                                     (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntries(workflowExecutionUuid: String,
                                      metadataKey: String,
                                      timeout: Duration)(implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntries(workflowExecutionUuid: String,
                                      callFullyQualifiedName: String,
                                      jobIndex: Option[Int],
                                      jobAttempt: Option[Int],
                                      timeout: Duration)(implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntries(workflowUuid: String,
                                      metadataKey: String,
                                      callFullyQualifiedName: String,
                                      jobIndex: Option[Int],
                                      jobAttempt: Option[Int],
                                      timeout: Duration)(implicit ec: ExecutionContext): Nothing = notImplemented()

    override def queryMetadataEntryWithKeyConstraints(workflowExecutionUuid: String,
                                             metadataKeysToFilterFor: List[String],
                                             metadataKeysToFilterAgainst: List[String],
                                             metadataJobQueryValue: MetadataJobQueryValue,
                                             timeout: Duration)
                                            (implicit ec: ExecutionContext): Nothing = notImplemented()

    override def summarizeIncreasing(summaryNameIncreasing: String,
                                     startMetadataKey: String,
                                     endMetadataKey: String,
                                     nameMetadataKey: String,
                                     statusMetadataKey: String,
                                     submissionMetadataKey: String,
                                     parentWorkflowIdKey: String,
                                     rootWorkflowIdKey: String,
                                     labelMetadataKey: String,
                                     limit: Int,
                                     buildUpdatedSummary:
                                     (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry])
                                       => WorkflowMetadataSummaryEntry)
                                    (implicit ec: ExecutionContext): Nothing = notImplemented()

    /**
      * Retrieves a window of summarizable metadata satisfying the specified criteria.
      *
      * @param buildUpdatedSummary Takes in the optional existing summary and the metadata, returns the new summary.
      * @return A `Future` with the maximum metadataEntryId summarized by the invocation of this method.
      */
    override def summarizeDecreasing(summaryNameDecreasing: String,
                                     summaryNameIncreasing: String,
                                     startMetadataKey: String,
                                     endMetadataKey: String,
                                     nameMetadataKey: String,
                                     statusMetadataKey: String,
                                     submissionMetadataKey: String,
                                     parentWorkflowIdKey: String,
                                     rootWorkflowIdKey: String,
                                     labelMetadataKey: String,
                                     limit: Int,
                                     buildUpdatedSummary:
                                     (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry])
                                       => WorkflowMetadataSummaryEntry)
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
  }
}

object WriteMetadataActorSpec {

  val IntermittentException = new Exception("Simulated Database Flakiness Exception") with NoStackTrace

  class BatchSizeCountingWriteMetadataActor(override val batchSize: Int,
                                            override val flushRate: FiniteDuration,
                                            override val serviceRegistryActor: ActorRef,
                                            override val threshold: Int) extends WriteMetadataActor(batchSize, flushRate, serviceRegistryActor, threshold) {

    var batchSizes: Vector[Int] = Vector.empty
    var failureCount: Int = 0

    override val recentArrivalThreshold = Some(100.millis)

    override def process(e: NonEmptyVector[MetadataWriteAction]) = {
      batchSizes = batchSizes :+ e.length
      val result = super.process(e)
      result.onComplete {
        case Success(_) => // Don't increment failure count
        case Failure(_) => failureCount += 1
      }

      result
    }
  }

}
