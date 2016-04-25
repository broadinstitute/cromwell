package cromwell.backend.impl.local

import akka.actor.ActorSystem
import akka.testkit.{ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import lenthall.exception.AggregatedException
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import spray.json.{JsObject, JsValue}
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlSource}

import scala.concurrent.duration._

class LocalInitializationActorSpec extends TestKit(ActorSystem("LocalInitializationActorSpec"))
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

  val HelloWorldTwoTasks =
    """
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello ${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |}
      |task ask {
      |  String message
      |  command {
      |    echo "Did you say ${message}!"
      |  }
      |  output {
      |    String question = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |  call ask {
      |  	input: message = hello.salutation
      |  }
      |}
    """.stripMargin

  val defaultBackendConfig = new BackendConfigurationDescriptor("config", ConfigFactory.load())

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

  private def getLocalBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(LocalInitializationActor.props(workflowDescriptor, calls, conf))
  }

  override def afterAll {
    system.shutdown()
  }

  "LocalInitializationActor" should {
    "return InitializationSuccess when there are no runtime attributes defined." in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid Docker entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid Docker entry based on input" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "\${addressee}" }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid Docker entry based on input from task 1" ignore {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorldTwoTasks, runtime = """runtime { docker: "\${message}" }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid Docker entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: 1 }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting docker runtime attribute to be a String"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting docker runtime attribute to a String'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid failOnStderr entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" failOnStderr: "false" }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid failOnStderr entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" failOnStderr: "yes" }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false''.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid continueOnReturnCode entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" failOnStderr: "false" continueOnReturnCode: 1}""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid continueOnReturnCode entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" failOnStderr: "yes" continueOnReturnCode: "value" }""")
        val backend = getLocalBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]''.")
            }
        }
      }
    }
  }
}
