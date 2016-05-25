package cromwell.backend.impl.local

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.io.BackendTestkitSpec
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s.Call

import scala.concurrent.duration._

class LocalInitializationActorSpec extends TestKit(ActorSystem("LocalInitializationActorSpec", ConfigFactory.parseString(
  """akka.loggers = ["akka.testkit.TestEventListener"]"""))) with BackendTestkitSpec with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender {

  val globalConfig = ConfigFactory.load()
  val backendConfig = globalConfig.getConfig("backend.providers.Local.config")
  val defaultBackendConfigDescriptor = new BackendConfigurationDescriptor(backendConfig, globalConfig)

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

  private def getLocalBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(LocalInitializationActor.props(workflowDescriptor, calls, conf))
  }

  override def afterAll {
    system.shutdown()
  }

  "LocalInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        EventFilter.warning(message = s"Key/s [memory] is/are not supported by LocalBackend. Unsupported attributes will not be part of jobs executions.", occurrences = 1) intercept {
          val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 }""")
          val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfigDescriptor)
          backend ! Initialize
        }
      }
    }
  }
}
