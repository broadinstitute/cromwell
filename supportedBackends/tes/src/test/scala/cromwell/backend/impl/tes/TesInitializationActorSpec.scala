package cromwell.backend.impl.tes

import java.util.UUID
import akka.actor.Props
import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.async.RuntimeAttributeValidationFailures
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.logging.LoggingTest._
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import spray.json.{JsNumber, JsObject, JsString}
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class TesInitializationActorSpec extends TestKitSuite
  with AnyWordSpecLike with Matchers with ImplicitSender {
  val Timeout: FiniteDuration = 10.second.dilated

  val HelloWorld: String =
    s"""
       |task hello {
       |  String addressee = "you"
       |  command {
       |    echo "Hello $${addressee}!"
       |  }
       |  output {
       |    String salutation = read_string(stdout())
       |  }
       |
       |  RUNTIME
       |}
       |
       |workflow wf_hello {
       |  call hello
       |}
    """.stripMargin

  val globalConfig: Config = ConfigFactory.parseString("")

  val backendConfigTemplate: String =
    """
      |// Base bucket for workflow executions
      |root = "cromwell-executions"
      |endpoint = "0.0.0.0"
      |
      |// Polling for completion backs-off gradually for slower-running jobs.
      |// This is the maximum polling interval (in seconds):
      |maximum-polling-interval = 600
      |
      |default-runtime-attributes {
      |    cpu: 1
      |    failOnStderr: false
      |    continueOnReturnCode: 0
      |    memory: "2 GB"
      |    disk: "2 GB"
      |    preemptible: false
      |    # The keys below have been commented out as they are optional runtime attributes.
      |    # dockerWorkingDir
      |    # docker
      |}
      |""".stripMargin


  private def getActorRef(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[CommandCallNode],
                          conf: BackendConfigurationDescriptor) = {
    val params = TesInitializationActorParams(workflowDescriptor, calls, new TesConfiguration(conf), emptyActor)
    val props = Props(new TesInitializationActor(params))
    system.actorOf(props, "TesInitializationActor" + UUID.randomUUID)
  }

  val backendConfig: Config = ConfigFactory.parseString(backendConfigTemplate)
  val conf: BackendConfigurationDescriptor = new BackendConfigurationDescriptor(backendConfig, globalConfig) {
    override private[backend] lazy val cromwellFileSystems = new CromwellFileSystems(globalConfig)
  }

  // TODO WOM: needs runtime attributes validation working again
  "TesInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld,
          runtime = """runtime { docker: "ubuntu/latest" test: true }""")
        val backend = getActorRef(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, conf)
        val eventPattern =
          "Key/s [test] is/are not supported by backend. Unsupported attributes will not be part of job executions."
        EventFilter.warning(pattern = escapePattern(eventPattern), occurrences = 1) intercept {
          backend ! Initialize
        }
        expectMsgPF() {
          case InitializationSuccess(_) => // Docker entry is present.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    def initializeActor(workflowOptions: WorkflowOptions): Unit = {
      val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld,
        runtime = """runtime { docker: "ubuntu/latest" }""",
        options = workflowOptions)
      val backend = getActorRef(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, conf)
      backend ! Initialize
    }

    def nonStringErrorMessage(key: String) = s"Workflow option $key must be a string"
    val bothRequiredErrorMessage = s"Workflow options ${TesWorkflowOptionKeys.WorkflowExecutionIdentity} and ${TesWorkflowOptionKeys.DataAccessIdentity} are both required if one is provided"

    "fail when WorkflowExecutionIdentity is not a string and DataAccessIdentity is missing" in {
      within(Timeout) {
        val workflowOptions = WorkflowOptions(
          JsObject(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsNumber(5)))
        )
        initializeActor(workflowOptions)
        expectMsgPF() {
          case InitializationSuccess(s) => fail(s"InitializationFailed was expected but got $s")
          case InitializationFailed(failure) =>
            val expectedMsg = nonStringErrorMessage(TesWorkflowOptionKeys.WorkflowExecutionIdentity)
            if (!(failure.getMessage.contains(expectedMsg) &&
                  failure.getMessage.contains(bothRequiredErrorMessage))) {
              fail(s"Exception message did not contain both '$expectedMsg' and '$bothRequiredErrorMessage'. Was '$failure'")
            }
        }
      }
    }

    "fail when WorkflowExecutionIdentity is a string but DataAccessIdentity is not a string" in {
      within(Timeout) {
        val workflowOptions = WorkflowOptions(JsObject(Map(
          TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsString("5"),
          TesWorkflowOptionKeys.DataAccessIdentity -> JsNumber(6)
        )))
        initializeActor(workflowOptions)
        expectMsgPF() {
          case InitializationSuccess(s) => fail(s"InitializationFailed was expected but got $s")
          case InitializationFailed(failure) =>
            val expectedMsg = nonStringErrorMessage(TesWorkflowOptionKeys.DataAccessIdentity)
            if (!failure.getMessage.contains(expectedMsg)) fail(s"Exception message did not contain '$expectedMsg'")
        }
      }
    }

    "successfully start when both WorkflowExecutionIdentity and DataAccessIdentity are strings" in {
      within(Timeout) {
        val workflowOptions = WorkflowOptions(JsObject(Map(
          TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsString("5"),
          TesWorkflowOptionKeys.DataAccessIdentity -> JsString("6")
        )))
        initializeActor(workflowOptions)
        expectMsgPF() {
          case InitializationSuccess(_) =>
          case InitializationFailed(f) => fail(s"InitializationSuccess was expected but got $f")
        }
      }
    }

    "return InitializationFailed when docker runtime attribute key is not present" in {
      within(Timeout) {
        val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
        val backend = getActorRef(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, conf)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: RuntimeAttributeValidationFailures =>
                if (!exception.getMessage.equals("Runtime validation failed:\nTask hello has an invalid runtime attribute docker = !! NOT FOUND !!"))
                  fail("Exception message is not equal to 'Runtime validation failed:\nTask hello has an invalid runtime attribute docker = !! NOT FOUND !!'.")
            }
        }
      }
    }
  }
}
