package cromwell.backend.impl.jes

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

class JesInitializationActorSpec extends TestKit(ActorSystem("JesInitializationActorSpec"))
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

  val defaultBackendConfig = new BackendConfigurationDescriptor("jes-config", ConfigFactory.load())

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
    "return InitializationFailed when there are no runtime attributes defined." in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
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

    "return InitializationSuccess when tries to validate a valid cpu entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" cpu: 1}""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid cpu entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" cpu: "value" }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting cpu runtime attribute to be an Integer"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting cpu runtime attribute to be an Integer'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid zones entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" zones: "us-central1-a" }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid zones entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" zones: 1 }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid array zones entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" zones: ["us-central1-a", "us-central1-b"] }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid array zones entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" zones: [2, 1] }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid preemptible entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" preemptible: 1}""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid preemptible entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" preemptible: "value" }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting preemptible runtime attribute to be an Integer"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting preemptible runtime attribute to be an Integer'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid bootDiskSizeGb entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" bootDiskSizeGb: 1}""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid bootDiskSizeGb entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" bootDiskSizeGb: "value" }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting bootDiskSizeGb runtime attribute to be an Integer"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting bootDiskSizeGb runtime attribute to be an Integer'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid disks entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" disks: "local-disk 10 SSD"}""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid disks entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" disks: 10 }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting disks runtime attribute to be a comma separated String or Array[String]"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting disks runtime attribute to be a comma separated String or Array[String]'.")
            }
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid array disks entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" disks: [10, 11] }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting disks runtime attribute to be a comma separated String or Array[String]"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting disks runtime attribute to be a comma separated String or Array[String]'.")
            }
        }
      }
    }

    "return InitializationSuccess when tries to validate a valid memory entry" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" memory: 1}""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationSuccess => //Entry is valid as expected.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when tries to validate an invalid memory entry" ignore {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu:latest" memory: "value" }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: AggregatedException =>
                if (!exception.exceptionContext.equals("Runtime attribute validation failed"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
                if (exception.throwables.size != 1)
                  fail("Number of errors coming from AggregatedException is not equals to one.")
                if (!exception.throwables.head.getMessage.contains("Expecting memory runtime attribute to be an Integer or String with format '8 GB'"))
                  fail("Message from error nr 1 in validation exception does not contains 'Expecting memory runtime attribute to be an Integer or String with format '8 GB''.")
            }
        }
      }
    }
  }
}
