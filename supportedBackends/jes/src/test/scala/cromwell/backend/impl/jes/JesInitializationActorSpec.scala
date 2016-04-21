package cromwell.backend.impl.jes

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import spray.json.{JsObject, JsValue}
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlSource}

import scala.concurrent.duration._

class JesInitializationActorSpec extends TestKit(ActorSystem("JesInitializationActorSpec", ConfigFactory.parseString(
  """akka.loggers = ["akka.testkit.TestEventListener"]"""))) with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender {
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

  private def getJesBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(JesInitializationActor.props(workflowDescriptor, calls, conf))
  }

  override def afterAll {
    system.shutdown()
  }

  "JesInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu/latest" test: true }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        EventFilter.warning(message = s"Key/s [test] is/are not supported by JesBackend. Unsupported attributes will not be part of jobs executions.", occurrences = 1) intercept {
          //Log message was intercepted.
        }
        expectMsgPF() {
          case InitializationSuccess => //Docker entry is present.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when docker runtime attribute key is not present" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: IllegalArgumentException =>
                if (!exception.getMessage.equals("docker mandatory runtime attribute is missing."))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
            }
        }
      }
    }
  }
}
