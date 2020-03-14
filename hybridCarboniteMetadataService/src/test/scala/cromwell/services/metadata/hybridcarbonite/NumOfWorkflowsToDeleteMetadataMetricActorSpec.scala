package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.testkit.{ImplicitSender, TestFSMRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.services.instrumentation.CromwellGauge
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.metadata.hybridcarbonite.NumOfWorkflowsToDeleteMetadataMetricActor.{CalculateNumOfWorkflowsToDeleteMetadataMetricValue, MetricCalculationInProgress, NumOfWorkflowsToDeleteMetadataMetricValue, WaitingForMetricCalculationRequestOrMetricValue}
import org.scalatest.{FlatSpecLike, Matchers}
import org.scalatest.concurrent.Eventually

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.concurrent.duration._
import scala.language.postfixOps

class NumOfWorkflowsToDeleteMetadataMetricActorSpec extends TestKitSuite with FlatSpecLike with ImplicitSender with Matchers with Eventually {

  private val defaultTimeout = 5 seconds

  it should "follow the proper state transition path when requested to calculate metric value and calculation succeeded" in {
    val serviceRegistryProbe = TestProbe()
    val metricActor = TestFSMRef(new NumOfWorkflowsToDeleteMetadataMetricActor(serviceRegistryProbe.ref))
    metricActor ! CalculateNumOfWorkflowsToDeleteMetadataMetricValue(OffsetDateTime.now())
    eventually {
      metricActor.stateName shouldBe MetricCalculationInProgress
    }
    eventually {
      metricActor.stateName shouldBe WaitingForMetricCalculationRequestOrMetricValue
    }
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        // real value should be zero in the empty test HSQLDB database
        actualValue shouldBe 0
    }
  }

  it should "return into WaitingForMetricCalculationRequestOrMetricValue state if DB query resulted in error" in {
    val dbFailurePromise = Promise[Int]()
    val metricActor = TestFSMRef(new NumOfWorkflowsToDeleteMetadataMetricActor(TestProbe().ref) {
      override def countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(archiveStatus: Option[String], thresholdTimestamp: OffsetDateTime)(implicit ec: ExecutionContext): Future[Int] =
        dbFailurePromise.future
    })
    metricActor ! CalculateNumOfWorkflowsToDeleteMetadataMetricValue(OffsetDateTime.now())
    eventually {
      metricActor.stateName shouldBe MetricCalculationInProgress
    }
    dbFailurePromise.failure(new RuntimeException("TEST: cannot query data from DB"))
    eventually {
      metricActor.stateName shouldBe WaitingForMetricCalculationRequestOrMetricValue
    }
  }

  it should "not query the database if received pre-calculated metric value and stay in the default state" in {
    val valFromDB = -1
    val precalculatedVal = -2L
    val serviceRegistryProbe = TestProbe()
    val dbResultPromise = Promise[Int]()
    val metricActor = TestFSMRef(new NumOfWorkflowsToDeleteMetadataMetricActor(serviceRegistryProbe.ref) {
      override def countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(archiveStatus: Option[String], thresholdTimestamp: OffsetDateTime)(implicit ec: ExecutionContext): Future[Int] = {
        dbResultPromise.success(valFromDB)
        dbResultPromise.future
      }
    })

    metricActor ! NumOfWorkflowsToDeleteMetadataMetricValue(precalculatedVal)
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        actualValue shouldBe precalculatedVal
    }

    //db query has not been called
    dbResultPromise.isCompleted shouldBe false
  }
}
