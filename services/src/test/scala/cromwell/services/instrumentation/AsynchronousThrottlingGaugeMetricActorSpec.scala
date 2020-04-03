package cromwell.services.instrumentation

import akka.testkit.{ImplicitSender, TestFSMRef, TestProbe}
import cats.data.NonEmptyList
import cromwell.core.TestKitSuite
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.instrumentation.AsynchronousThrottlingGaugeMetricActor.{CalculateMetricValue, MetricCalculationInProgress, MetricValue, WaitingForMetricCalculationRequestOrMetricValue}
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import org.scalatest.{FlatSpecLike, Matchers}
import org.scalatest.concurrent.Eventually

import scala.concurrent.Promise
import scala.concurrent.duration._
import scala.language.postfixOps

class AsynchronousThrottlingGaugeMetricActorSpec extends TestKitSuite with FlatSpecLike with ImplicitSender with Matchers with Eventually {

  private val defaultTimeout = 5 seconds

  it should "follow the proper state transition path when requested to calculate metric value and calculation succeeded" in {
    val calculatedVal = -1
    val serviceRegistryProbe = TestProbe()
    val calculatedValPromise = Promise[Int]()
    val metricActor = TestFSMRef {
      new AsynchronousThrottlingGaugeMetricActor(NonEmptyList.of("metric"), InstrumentationPrefixes.ServicesPrefix, serviceRegistryProbe.ref)
    }
    metricActor ! CalculateMetricValue(_ => calculatedValPromise.future)
    eventually {
      metricActor.stateName shouldBe MetricCalculationInProgress
    }
    calculatedValPromise.success(calculatedVal)
    eventually {
      metricActor.stateName shouldBe WaitingForMetricCalculationRequestOrMetricValue
    }
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        actualValue shouldBe calculatedVal
    }
  }

  it should "return into WaitingForMetricCalculationRequestOrMetricValue state if metric calculation resulted in error" in {
    val dbFailurePromise = Promise[Int]()
    val metricActor = TestFSMRef {
      new AsynchronousThrottlingGaugeMetricActor(NonEmptyList.of("metric"), InstrumentationPrefixes.ServicesPrefix, TestProbe().ref)
    }
    metricActor ! CalculateMetricValue(_ => dbFailurePromise.future)
    eventually {
      metricActor.stateName shouldBe MetricCalculationInProgress
    }
    dbFailurePromise.failure(new RuntimeException("TEST: cannot calculate metric value"))
    eventually {
      metricActor.stateName shouldBe WaitingForMetricCalculationRequestOrMetricValue
    }
  }

  it should "submit metric value to instrumentation and stay in WaitingForMetricCalculationRequestOrMetricValue state if received pre-calculated metric value" in {
    val calculatedVal = -1L
    val serviceRegistryProbe = TestProbe()
    val metricActor = TestFSMRef {
      new AsynchronousThrottlingGaugeMetricActor(NonEmptyList.of("metric"), InstrumentationPrefixes.ServicesPrefix, serviceRegistryProbe.ref)
    }
    metricActor ! MetricValue(calculatedVal)
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        actualValue shouldBe calculatedVal
    }
  }

  it should "successfully complete ongoing metric value calculation even after being interrupted by MetricValue message" in {
    val calculatedVal = -1
    val precalculatedVal = -2L
    val serviceRegistryProbe = TestProbe()
    val calculatedValPromise = Promise[Int]()
    val metricActor = TestFSMRef {
      new AsynchronousThrottlingGaugeMetricActor(NonEmptyList.of("metric"), InstrumentationPrefixes.ServicesPrefix, serviceRegistryProbe.ref)
    }
    metricActor ! CalculateMetricValue(_ => calculatedValPromise.future)
    eventually {
      metricActor.stateName shouldBe MetricCalculationInProgress
    }

    // interrupted by precalculated MetricValue
    metricActor ! MetricValue(precalculatedVal)
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        actualValue shouldBe precalculatedVal
    }
    metricActor.stateName shouldBe MetricCalculationInProgress

    // completing ongoing calculation
    calculatedValPromise.success(calculatedVal)
    eventually {
      metricActor.stateName shouldBe WaitingForMetricCalculationRequestOrMetricValue
    }
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        actualValue shouldBe calculatedVal
    }
  }

  it should "successfully complete ongoing metric value calculation and ignore another calculation requests while in MetricCalculationInProgress state" in {
    val serviceRegistryProbe = TestProbe()
    val metricActor = TestFSMRef {
      new AsynchronousThrottlingGaugeMetricActor(NonEmptyList.of("metric"), InstrumentationPrefixes.ServicesPrefix, serviceRegistryProbe.ref)
    }

    val calculatedVal = -1
    val calculatedValPromise = Promise[Int]()
    val calculatedValInterruptor = -2
    val calculatedValInterruptorPromise = Promise[Int]()

    metricActor ! CalculateMetricValue(_ => calculatedValPromise.future)
    eventually {
      metricActor.stateName shouldBe MetricCalculationInProgress
    }

    metricActor ! CalculateMetricValue(_ => calculatedValInterruptorPromise.future)

    // completing ongoing calculation
    calculatedValPromise.success(calculatedVal)
    // complete interrupttor calculation
    calculatedValInterruptorPromise.success(calculatedValInterruptor)
    eventually {
      metricActor.stateName shouldBe WaitingForMetricCalculationRequestOrMetricValue
    }
    serviceRegistryProbe.expectMsgPF(defaultTimeout){
      case InstrumentationServiceMessage(CromwellGauge(_, actualValue)) =>
        // should be calculatedVal, not calculatedValInterruptor
        actualValue shouldBe calculatedVal
    }
    serviceRegistryProbe.expectNoMessage(defaultTimeout)
  }
}
