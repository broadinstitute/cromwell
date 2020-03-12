package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.Props
import akka.testkit.{ImplicitSender, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor.MetricActorFreed
import cromwell.services.metadata.hybridcarbonite.WorkflowsToDeleteMetadataMetricHelperActor.{CalculateMetricAndSend, SendPrecalculatedMetric}
import org.scalatest.FlatSpecLike

import scala.concurrent.{ExecutionContext, Future}

class WorkflowsToDeleteMetadataMetricHelperActorSpec extends TestKitSuite with FlatSpecLike with ImplicitSender {

  it should "send MetricActorFreed to DeleteMetadataActor after restart" in {
    val deleteMetadataActor = TestProbe()
    val metricActor = system.actorOf(Props(new WorkflowsToDeleteMetadataMetricHelperActor(deleteMetadataActor.ref, TestProbe().ref) {
      override def receive: Receive = {
        case SendPrecalculatedMetric(_) =>
          throw new RuntimeException("TEST: runtime exception")
      }
    }))
    metricActor ! SendPrecalculatedMetric(-1L)
    deleteMetadataActor.expectMsg(MetricActorFreed)
  }

  it should "send MetricActorFreed to DeleteMetadataActor after sending metric completed with either Success or Failure" in {
    val deleteMetadataActor = TestProbe()

    val metricActorSuccessful = system.actorOf(Props(new WorkflowsToDeleteMetadataMetricHelperActor(deleteMetadataActor.ref, TestProbe().ref) {
      override def countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(archiveStatus: Option[String], thresholdTimestamp: OffsetDateTime)(implicit ec: ExecutionContext): Future[Int] = {
        Future.successful(0)
      }
    }))
    val metricActorFailing = system.actorOf(Props(new WorkflowsToDeleteMetadataMetricHelperActor(deleteMetadataActor.ref, TestProbe().ref) {
      override def countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(archiveStatus: Option[String], thresholdTimestamp: OffsetDateTime)(implicit ec: ExecutionContext): Future[Int] = {
        Future.failed(new RuntimeException("TEST: cannot calculate metric"))
      }
    }))

    metricActorSuccessful ! SendPrecalculatedMetric(-1L)
    deleteMetadataActor.expectMsg(MetricActorFreed)

    metricActorSuccessful ! CalculateMetricAndSend(OffsetDateTime.now())
    deleteMetadataActor.expectMsg(MetricActorFreed)

    metricActorFailing ! CalculateMetricAndSend(OffsetDateTime.now())
    deleteMetadataActor.expectMsg(MetricActorFreed)
  }

}
