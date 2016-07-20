package cromwell.services

import akka.actor.ActorSystem
import akka.testkit.TestKit
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import akka.testkit._

import scala.concurrent.duration._

object ServicesSpec {
  val ConfigText =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "INFO"

      |  test {
      |    # Some of our tests fire off a message, then expect a particular event message within 3s (the default).
      |    # Especially on CI, the metadata test does not seem to be returning in time. So, overriding the timeouts
      |    # with slightly higher values. Alternatively, could also adjust the akka.test.timefactor only in CI.
      |    filter-leeway = 5s
      |    single-expect-default = 5s
      |    default-timeout = 10s
      |  }
      |}
    """.stripMargin

  val config = ConfigFactory.parseString(ServicesSpec.ConfigText)
}

abstract class ServicesSpec(serviceName: String) extends TestKit(ActorSystem(s"${serviceName}ServiceActorSpec", ServicesSpec.config))
  with Matchers with WordSpecLike with BeforeAndAfterAll with ScalaFutures {

  implicit val timeout = Timeout(10.seconds.dilated)
  implicit val ec = system.dispatcher
  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  override protected def afterAll() = {
    TestKit.shutdownActorSystem(system)
  }
}
