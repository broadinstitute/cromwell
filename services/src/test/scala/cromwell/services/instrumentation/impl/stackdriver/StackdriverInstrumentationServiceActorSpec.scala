package cromwell.services.instrumentation.impl.stackdriver

import akka.testkit.{TestActorRef, TestProbe}
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.instrumentation._
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._

class StackdriverInstrumentationServiceActorSpec extends TestKitSuite with FlatSpecLike with Matchers with Eventually {
  behavior of "StackdriverInstrumentationServiceActor"

  val MaxWaitTime = 2.minutes
  implicit val pc: PatienceConfig = PatienceConfig(MaxWaitTime)

  val globalConfig = ConfigFactory.parseString(
    s"""
       |google {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = "application-default"
       |      scheme = "application_default"
       |    }
       |  ]
       |}
      """.stripMargin
  )

  val registryProbe = TestProbe().ref
  val resourceLabels = Map("project_id" -> "my-project")

  val testBucket = CromwellBucket(List("test_prefix"), NonEmptyList.of("test", "metric", "bucket"))
  val testGaugeBucket = CromwellBucket(List("test_prefix"), NonEmptyList.of("test", "gauge", "metric", "bucket"))


  it should "correctly receive the metrics with resource labels" in {
    val stackdriverConfig = ConfigFactory.parseString(
      """
        |auth = "application-default"
        |google-project = "my-project"
        |flush-rate = 1 minute
      """.stripMargin
    )

    val rawMetricList = List(
      CromwellIncrement(testBucket),
      CromwellCount(testBucket, 80),
      CromwellGauge(testGaugeBucket, 20),
      CromwellTiming(testBucket.expand("timing"), 5.seconds),
      CromwellIncrement(testBucket),
      CromwellCount(testBucket, 80),
      CromwellGauge(testGaugeBucket, 60),
      CromwellTiming(testBucket.expand("timing"), 10.seconds)
    )

    val expectedCumulativeMetrics = ("custom.googleapis.com/cromwell/test_prefix/test/metric/bucket", 162.0)
    val expectedGaugeMetrics = ("custom.googleapis.com/cromwell/test_prefix/test/gauge/metric/bucket", 40.0)
    val expectedTimingMetrics = ("custom.googleapis.com/cromwell/test_prefix/test/metric/bucket/timing", 7500.0)

    val stackdriverActor = TestActorRef(new TestStackdriverInstrumentationServiceActor(stackdriverConfig, globalConfig, registryProbe))

    rawMetricList foreach (metric => stackdriverActor ! InstrumentationServiceMessage(metric))

    eventually {
      stackdriverActor.underlyingActor.metricsReceived should have length 3

      stackdriverActor.underlyingActor.metricsReceived.map(m => (m.metricPath, m.metricValue)) should contain (expectedCumulativeMetrics)
      stackdriverActor.underlyingActor.metricsReceived.map(m => (m.metricPath, m.metricValue)) should contain (expectedGaugeMetrics)
      stackdriverActor.underlyingActor.metricsReceived.map(m => (m.metricPath, m.metricValue)) should contain (expectedTimingMetrics)

      stackdriverActor.underlyingActor.metricsReceived.map(_.resourceLabels) should contain (resourceLabels)
    }
  }


  it should "correctly receive metrics with metric labels" in {
    val stackdriverConfig = ConfigFactory.parseString(
      """
        |auth = "application-default"
        |google-project = "my-project"
        |flush-rate = 1 minute
        |cromwell-instance-role = "backend"
        |cromwell-perf-test-case = "perf-test-1"
      """.stripMargin
    )

    val rawMetricList = List(
      CromwellIncrement(testBucket),
      CromwellCount(testBucket, 49)
    )

    val expectedCumulativeMetrics = ("custom.googleapis.com/cromwell/test_prefix/test/metric/bucket", 50.0)
    val metricLabels = Map("cromwell_instance_role" -> "backend", "cromwell_perf_test_case" -> "perf-test-1")

    val stackdriverActor = TestActorRef(new TestStackdriverInstrumentationServiceActor(stackdriverConfig, globalConfig, registryProbe))

    rawMetricList foreach (metric => stackdriverActor ! InstrumentationServiceMessage(metric))

    eventually {
      stackdriverActor.underlyingActor.metricsReceived should have length 1
      stackdriverActor.underlyingActor.metricsReceived.map(m => (m.metricPath, m.metricValue)) should contain (expectedCumulativeMetrics)
      stackdriverActor.underlyingActor.metricsReceived.map(_.resourceLabels) should contain (resourceLabels)
      stackdriverActor.underlyingActor.metricsReceived.map(_.metricLabels) should contain (metricLabels)
    }
  }
}
