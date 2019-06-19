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
    """
      |hostname = "localhost"
      |port = 8125
      |prefix = "prefix_value" # can be used to prefix all metrics with an api key for example
      |flush-rate = 100 ms # rate at which aggregated metrics will be sent to statsd
    """.stripMargin
  )

  val cromwellInstanceConfig = ConfigFactory.parseString("""system.cromwell_id = "cromwell-101"""")

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

  // Fuzzy packets are metrics which we expect to see but for which cannot determine the exact values for asynchronous reasons
  case class StatsDTestBit(description: String,
                           metric: CromwellMetric,
                           expectedExactPackets: Set[String],
                           expectedFuzzyPackets: Set[String] = Set.empty,
                           prependCromwellInstancePrefix: Boolean = false)
  
  // Note: The current StatsD implementation sends everything as StatsD gauges so we expect all packets to be "...|g"
  List(
    StatsDTestBit("increment counters", CromwellIncrement(testBucket),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.samples:1|g"),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.m1_rate")
    ),
    StatsDTestBit("add count", CromwellCount(testBucket, 80),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.samples:81|g"),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.m1_rate")
    ),
    StatsDTestBit("set gauges", CromwellGauge(testGaugeBucket, 89),
      Set("prefix_value.cromwell.test_prefix.test.gauge.metric.bucket:89|g")
    ),
    StatsDTestBit("set timings", CromwellTiming(testBucket.expand("timing"), 5.seconds),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.timing.stddev:0.00|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.samples:1|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.p95:5000.00|g",
        "prefix_value.cromwell.test_prefix.test.metric.bucket.timing.mean:5000.00|g"
      ),
      Set("prefix_value.cromwell.test_prefix.test.metric.bucket.timing.m1_rate")
    ),
    StatsDTestBit("increment counters with cromwell id prefix", CromwellIncrement(testBucket),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.samples:1|g"),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.m1_rate"),
      prependCromwellInstancePrefix = true
    ),
    StatsDTestBit("add count with cromwell id prefix", CromwellCount(testBucket, 80),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.samples:81|g"),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.m1_rate"),
      prependCromwellInstancePrefix = true
    ),
    StatsDTestBit("set gauges with cromwell id prefix", CromwellGauge(testGaugeBucket, 89),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.gauge.metric.bucket:89|g"),
      prependCromwellInstancePrefix = true
    ),
    StatsDTestBit("set timings with cromwell id prefix", CromwellTiming(testBucket.expand("timing"), 5.seconds),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.timing.stddev:0.00|g",
        "prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.timing.samples:1|g",
        "prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.timing.p95:5000.00|g",
        "prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.timing.mean:5000.00|g"
      ),
      Set("prefix_value.cromwell.cromwell-101.test_prefix.test.metric.bucket.timing.m1_rate"),
      prependCromwellInstancePrefix = true
    )
  ) foreach {
    case StatsDTestBit(description, metric, expectedExactPackets, expectedFuzzyPackets, prependCromwellInstancePrefix) =>
      it should description in {

        val globalConfig = if (prependCromwellInstancePrefix) cromwellInstanceConfig else ConfigFactory.empty()

        val instrumentationActor = TestActorRef(new StatsDInstrumentationServiceActor(config, globalConfig, registryProbe))
        instrumentationActor ! InstrumentationServiceMessage(metric)
        val received = udpProbe.receiveWhile(patience) {
          case Udp.Received(data, _) => data.utf8String
        }

        expectedExactPackets foreach { packet => if (!received.contains(packet)) {
          val prefix = packet.split(":").head
          received.find(_.startsWith(prefix)) match {
            case Some(sharedPrefix) => fail(s"Missing packet: $packet, but found: $sharedPrefix. Should this be a fuzzy packet?")
            case None => fail(s"Missing packet: $packet, and no packets received with prefix $prefix")
          }
        }}
        expectedFuzzyPackets foreach { packet => if (!received.exists(_.contains(packet))) fail(s"Missing fuzzy packet: $packet") }
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
