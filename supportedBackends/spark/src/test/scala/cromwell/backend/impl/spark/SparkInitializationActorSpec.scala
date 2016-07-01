package cromwell.backend.impl.spark

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendWorkflowDescriptor, BackendConfigurationDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import spray.json.{JsValue, JsObject}
import wdl4s._
import wdl4s.values.WdlValue
import scala.concurrent.duration._

class SparkInitializationActorSpec  extends TestKit(ActorSystem("SparkInitializationActorSpec", ConfigFactory.parseString(
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

  override def afterAll {
    system.shutdown()
  }

  private def getHtCondorBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(SparkInitializationActor.props(workflowDescriptor, calls, conf))
  }

  "SparkInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        EventFilter.warning(message = s"Key/s [memory] is/are not supported by SparkBackend. Unsupported attributes will not be part of jobs executions.", occurrences = 1) intercept {
          val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 }""")
          val backend = getHtCondorBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
          backend ! Initialize
        }
      }
    }
  }
}
