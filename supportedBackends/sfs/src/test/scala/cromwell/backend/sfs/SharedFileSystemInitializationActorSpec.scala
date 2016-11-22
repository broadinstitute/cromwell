package cromwell.backend.sfs

import akka.actor.Props
import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.TestKitSuite
import cromwell.core.logging.LoggingTest._
import org.scalatest.{Matchers, WordSpecLike}
import wdl4s.TaskCall

import scala.concurrent.duration._

class SharedFileSystemInitializationActorSpec extends TestKitSuite("SharedFileSystemInitializationActorSpec")
  with WordSpecLike with Matchers with ImplicitSender {
  val Timeout = 5.second.dilated

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

  private def getActorRef(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                          conf: BackendConfigurationDescriptor) = {
    val params = SharedFileSystemInitializationActorParams(emptyActor, workflowDescriptor, conf, calls, List.empty)
    val props = Props(new SharedFileSystemInitializationActor(params))
    system.actorOf(props, "SharedFileSystemInitializationActor")
  }

  "SharedFileSystemInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { unsupported: 1 }""")
        val conf = emptyBackendConfig
        val backend = getActorRef(workflowDescriptor, workflowDescriptor.workflow.taskCalls, conf)
        val pattern = "Key/s [unsupported] is/are not supported by backend. " +
          "Unsupported attributes will not be part of jobs executions."
        EventFilter.warning(pattern = escapePattern(pattern), occurrences = 1) intercept {
          backend ! Initialize
        }
      }
    }
  }
}
