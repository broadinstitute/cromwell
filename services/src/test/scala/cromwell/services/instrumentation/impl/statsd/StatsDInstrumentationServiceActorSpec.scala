package cromwell.services.instrumentation.impl.statsd

import java.net.InetSocketAddress

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.io.{IO, Udp}
import akka.testkit.{TestActorRef, TestProbe}
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.instrumentation._
import org.scalatest.{BeforeAndAfterAll, FlatSpecLike, Matchers}

import scala.concurrent.duration._
class StatsDInstrumentationServiceActorSpec extends TestKitSuite with FlatSpecLike with BeforeAndAfterAll with Matchers {
  behavior of "StatsDInstrumentationServiceActor"

  val config = ConfigFactory.parseString(
    """statsd {
      |hostname = "localhost"
      |port = 8125
      |prefix = "prefix_value" # can be used to prefix all metrics with an api key for example
      |flush-rate = 100 ms # rate at which aggregated metrics will be sent to statsd
      |}
    """.stripMargin
  )

  val registryProbe = TestProbe().ref
  val udpProbe = TestProbe()
  val patience = 1.second
  val testBucket = CromwellBucket(List("test_prefix"), NonEmptyList.of("test", "metric", "bucket"))
  val testGaugeBucket = CromwellBucket(List("test_prefix"), NonEmptyList.of("test", "gauge", "metric", "bucket"))
  
  var udpListenerActor: ActorRef = _
  
  override def beforeAll(): Unit = {
    // Start an actor listening to the UDP port and forwarding messages to the udpProbe
    udpListenerActor = system.actorOf(Props(new UDPListenerActor(new InetSocketAddress("localhost", 8125), udpProbe.ref)))
    // Give a sec to the actor to open an UDP socket
    Thread.sleep(1.second.toMillis)
    super.beforeAll()
  }

  override def afterAll(): Unit = {
    udpListenerActor ! Udp.Unbind
    super.afterAll()
  }

  case class StatsDTestBit(description: String, metric: CromwellMetric, expectedPackets: Set[String])
  
  // Note: The current StatsD implementation sends everything as StatsD gauges so we expect all packets to be "...|g"
  List(
    StatsDTestBit("increment counters", CromwellIncrement(testBucket),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.samples:1|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.m1_rate:0.00|g"
      )
    ),
    StatsDTestBit("add count", CromwellCount(testBucket, 80),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.samples:81|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.m1_rate:0.00|g"
      )
    ),
    StatsDTestBit("set gauges", CromwellGauge(testGaugeBucket, 89),
      Set("prefix_value.cromwell.test_prefix.test.gauge.metric.bucket:89|g")
    ),
    StatsDTestBit("set timings", CromwellTiming(testBucket.expand("timing"), 5.seconds),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.timing.stddev:0.00|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.samples:1|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.p95:5000.00|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.mean:5000.00|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.m1_rate:0.00|g"
      )
    )
  ) foreach {
    case StatsDTestBit(description, metric, expectedPackets) =>
      it should description in {
        val instrumentationActor = TestActorRef(new StatsDInstrumentationServiceActor(config, ConfigFactory.load(), registryProbe))
        instrumentationActor ! InstrumentationServiceMessage(metric)
        val received = udpProbe.receiveWhile(patience) {
          case Udp.Received(data, _) => data.utf8String
        }

        expectedPackets foreach { packet => if (!received.contains(packet)) fail(s"Missing packet: $packet") }
      }
  }
  
  private class UDPListenerActor(remote: InetSocketAddress, sendTo: ActorRef) extends Actor with ActorLogging {
    implicit val system =  context.system
    IO(Udp) ! Udp.Bind(sendTo, remote)

    def receive = {
      case Udp.Bound(_) => context.become(ready(sender()))
    }

    def ready(socket: ActorRef): Receive = {
      case Udp.Unbind  => socket ! Udp.Unbind
      case Udp.Unbound => context.stop(self)
      case other => log.error(s"received unexpected message: $other")
    }
  }
}
