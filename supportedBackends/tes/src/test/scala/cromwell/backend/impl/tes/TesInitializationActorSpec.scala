package cromwell.backend.impl.tes

import java.util.UUID

import akka.actor.Props
import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.async.RuntimeAttributeValidationFailures
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.Tags.PostWomTest
import cromwell.core.TestKitSuite
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.logging.LoggingTest._
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class TesInitializationActorSpec extends TestKitSuite("TesInitializationActorSpec")
  with AnyWordSpecLike with Matchers with ImplicitSender {
  val Timeout = 10.second.dilated

  val HelloWorld =
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
  val conf = new BackendConfigurationDescriptor(backendConfig, globalConfig) {
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

    "return InitializationFailed when docker runtime attribute key is not present" taggedAs PostWomTest ignore {
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

