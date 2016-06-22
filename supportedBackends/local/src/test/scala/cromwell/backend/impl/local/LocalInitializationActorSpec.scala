package cromwell.backend.impl.local

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.ConfigResourceString.usingAsConfigFile
import cromwell.backend.io.BackendTestkitSpec
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, ConfigResourceString}
import cromwell.core.logging.LoggingTest._
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s.Call

import scala.concurrent.duration._

class LocalInitializationActorSpec extends TestKit(ActorSystem("LocalInitializationActorSpec", ConfigFactory.parseString(
  // TODO: PBE: 5s leeway copy of CromwellTestkitSpec. Refactor to D.R.Y. this code, and rename Testkit to TestKit
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
    |""".stripMargin)))
  with BackendTestkitSpec with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender with HasLocalBackendConfig {

  val defaultBackendConfigDescriptor = new BackendConfigurationDescriptor(localBackendConfig, globalConfig = null)

  val Timeout = 5.seconds.dilated

  val HelloWorld =
    """
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello ${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin

  private def getLocalBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(LocalInitializationActor.props(workflowDescriptor, calls, conf))
  }

  override def afterAll {
    system.shutdown()
  }

  "LocalInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      usingAsConfigFile("services-application.conf") {
        within(Timeout) {
          EventFilter.warning(pattern = escapePattern(s"Key/s [memory] is/are not supported by LocalBackend. Unsupported attributes will not be part of jobs executions."), occurrences = 1) intercept {
            val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 }""")
            val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfigDescriptor)
            backend ! Initialize
          }
        }
      }
    }
  }
}
