package cromwell.backend.sfs

import akka.actor.Props
import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.standard.DefaultInitializationActorParams
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, TestConfig}
import cromwell.core.TestKitSuite
import cromwell.core.logging.LoggingTest._
import org.scalatest.{Matchers, WordSpecLike}
import wdl4s.TaskCall

import scala.concurrent.duration._

class SharedFileSystemInitializationActorSpec extends TestKitSuite("SharedFileSystemInitializationActorSpec")
  with WordSpecLike with Matchers with ImplicitSender {
  val Timeout: FiniteDuration = 5.second.dilated

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

  private def getActorRef(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                          conf: BackendConfigurationDescriptor) = {
    val params = DefaultInitializationActorParams(workflowDescriptor, emptyActor, calls, emptyActor, conf)
    val props = Props(new SharedFileSystemInitializationActor(params))
    system.actorOf(props, "SharedFileSystemInitializationActor")
  }

  "SharedFileSystemInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { unsupported: 1 }""")
        val conf = BackendConfigurationDescriptor(TestConfig.mockBackendRuntimeConfig, ConfigFactory.empty())
        val backend = getActorRef(workflowDescriptor, workflowDescriptor.workflow.taskCalls, conf)
        val pattern = "Key/s [unsupported] is/are not supported by backend. " +
          "Unsupported attributes will not be part of job executions."
        EventFilter.warning(pattern = escapePattern(pattern), occurrences = 1) intercept {
          backend ! Initialize
        }
      }
    }
  }
}
