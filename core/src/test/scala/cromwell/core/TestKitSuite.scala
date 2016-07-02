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
      |    slow-actor-dispatcher {
      |      type = Dispatcher
      |      executor = "fork-join-executor"
      |    }
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
