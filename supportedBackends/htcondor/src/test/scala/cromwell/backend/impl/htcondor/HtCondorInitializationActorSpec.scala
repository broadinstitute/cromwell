package cromwell.backend.impl.htcondor

import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendConfigurationDescriptor, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.TestKitSuite
import org.scalatest.{Matchers, WordSpecLike}
import wdl4s.Call

import scala.concurrent.duration._

class HtCondorInitializationActorSpec extends TestKitSuite("HtCondorInitializationActorSpec") with WordSpecLike
  with Matchers with ImplicitSender {
  val Timeout = 5.second.dilated

  import BackendSpec._

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

  private def getHtCondorBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(HtCondorInitializationActor.props(workflowDescriptor, calls, conf))
  }

  "HtCondorInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        EventFilter.warning(message = s"Key/s [proc] is/are not supported by HtCondorBackend. Unsupported attributes will not be part of jobs executions.", occurrences = 1) intercept {
          val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { proc: 1 }""")
          val backend = getHtCondorBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls,
            emptyBackendConfig)
          backend ! Initialize
        }
      }
    }
  }
}
