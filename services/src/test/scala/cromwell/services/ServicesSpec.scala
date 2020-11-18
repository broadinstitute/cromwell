package cromwell.services

import akka.testkit._
import akka.util.Timeout
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.matchers.should.Matchers
import org.scalatest.time.{Millis, Minute, Span}
import org.scalatest.wordspec.AnyWordSpecLike

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

object ServicesSpec {
  val configString: String =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "DEBUG"
      |  actor {
      |    debug {
      |       receive = on
      |    }
      |    guardian-supervisor-strategy = "akka.actor.DefaultSupervisorStrategy"
      |  }
      |  dispatchers {
      |    # A dispatcher for actors performing blocking io operations
      |    # Prevents the whole system from being slowed down when waiting for responses from external resources for instance
      |    io-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |      # Using the forkjoin defaults, this can be tuned if we wish
      |    }
      |
      |    # A dispatcher for actors handling API operations
      |    # Keeps the API responsive regardless of the load of workflows being run
      |    api-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |
      |    # A dispatcher for engine actors
      |    # Because backends behavior is unpredictable (potentially blocking, slow) the engine runs
      |    # on its own dispatcher to prevent backends from affecting its performance.
      |    engine-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |
      |    # A dispatcher used by supported backend actors
      |    backend-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
      |
      |    # Note that without further configuration, backend actors run on the default dispatcher
      |  }
      |  test {
      |    # Some of our tests fire off a message, then expect a particular event message within 3s (the default).
      |    # Especially on CI, the metadata test does not seem to be returning in time. So, overriding the timeouts
      |    # with slightly higher values. Alternatively, could also adjust the akka.test.timefactor only in CI.
      |    filter-leeway = 5s
      |    single-expect-default = 5s
      |    default-timeout = 10s
      |  }
      |}
      |""".stripMargin

  val config: Config = ConfigFactory.parseString(ServicesSpec.configString)
}

abstract class ServicesSpec extends TestKitSuite
  with Matchers with AnyWordSpecLike with ScalaFutures {

  override protected lazy val actorSystemConfig: Config = ServicesSpec.config
  implicit val timeout: Timeout = Timeout(20.seconds.dilated)
  implicit val ec: ExecutionContext = system.dispatcher
  implicit val defaultPatience: PatienceConfig = PatienceConfig(timeout = Span(1, Minute), interval = Span(100, Millis))
}
