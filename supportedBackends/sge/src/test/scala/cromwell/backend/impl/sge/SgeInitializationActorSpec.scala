package cromwell.backend.impl.sge

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import cromwell.core.logging.LoggingTest._
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import spray.json.{JsObject, JsValue}
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlSource}

import scala.concurrent.duration._

class SgeInitializationActorSpec extends TestKit(ActorSystem("SgeInitializationActorSpec", ConfigFactory.parseString(
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
  with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender {
  val Timeout = 5.second.dilated

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

  val defaultBackendConfig = new BackendConfigurationDescriptor(ConfigFactory.parseString("{}"), ConfigFactory.load())

  private def buildWorkflowDescriptor(wdl: WdlSource,
                                      inputs: Map[String, WdlValue] = Map.empty,
                                      options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                                      runtime: String = "") = {
    new BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime)),
      inputs,
      options
    )
  }

  private def getSgeBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(SgeInitializationActor.props(workflowDescriptor, calls, conf))
  }

  override def afterAll {
    system.shutdown()
  }

  "SgeInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 }""")
        val backend = getSgeBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        EventFilter.warning(pattern = escapePattern(s"Key/s [memory] is/are not supported by SgeBackend. Unsupported attributes will not be part of jobs executions."), occurrences = 1) intercept {
          backend ! Initialize
        }
      }
    }
  }
}
