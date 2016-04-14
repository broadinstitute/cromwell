import akka.actor.ActorSystem
import akka.testkit.{ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, Initialize}
import cromwell.backend.impl.htcondor.HtCondorInitializationActor
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import lenthall.exception.AggregatedException
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import spray.json.{JsObject, JsValue}
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlSource}

import scala.concurrent.duration._

class HtCondorInitializationActorSpec extends TestKit(ActorSystem("HtCondorInitializationActorSpec"))
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

  val defaultBackendConfig = new BackendConfigurationDescriptor("htcondor-config", ConfigFactory.load())

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
    system.actorOf(HtCondorInitializationActor.props(workflowDescriptor, calls, conf))
  }

  "HtCondorInitializationActor" should {
    "return InitializationFailed when there are no runtime attributes defined." in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
        val backend = getHtCondorBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Failed to get Docker mandatory key from runtime attributes"))
                  fail("Message from error nr 1 in validation exception does not contains 'Failed to get Docker mandatory key from runtime attributes'.")
            }
        }
      }
    }
  }
}
