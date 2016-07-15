package cromwell.core

import java.util.UUID

import akka.actor.ActorSystem
import akka.testkit.TestKit
import com.typesafe.config.{Config, ConfigFactory}
import org.scalatest.{BeforeAndAfterAll, Suite}

/**
  * A mix of Akka TestKit with ScalaTest mixed in to clean up the actor system.
  *
  * @param actorSystemName The name of the actor system.
  * @param actorSystemConfig The config for the actor system.
  */
abstract class TestKitSuite(actorSystemName: String = TestKitSuite.randomName,
                            actorSystemConfig: Config = TestKitSuite.config)
  extends TestKit(ActorSystem(actorSystemName, actorSystemConfig)) with Suite with BeforeAndAfterAll {

  override protected def afterAll() = {
    system.shutdown()
    super.afterAll()
  }
}

object TestKitSuite {
  val configString =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "INFO"
      |  actor {
      |    debug {
      |       receive = on
      |    }
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
      |    # Because backends behaviour is unpredictable (potentially blocking, slow) the engine runs
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

  val config = ConfigFactory.parseString(configString)

  def randomName = s"TestSystem-${UUID.randomUUID}"
}
