package cromwell.backend.impl.sge

import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendConfigurationDescriptor, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.TestKitSuite
import cromwell.core.logging.LoggingTest._
import org.scalatest.{Matchers, WordSpecLike}
import wdl4s.Call

import scala.concurrent.duration._

class SgeInitializationActorSpec extends TestKitSuite("SgeInitializationActorSpec") with WordSpecLike with Matchers
  with ImplicitSender {
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

  private def getSgeBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(SgeInitializationActor.props(workflowDescriptor, calls, conf, emptyActor))
  }

  "SgeInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 }""")
        val backend = getSgeBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls,
          emptyBackendConfig)
        EventFilter.warning(pattern = escapePattern(s"Key/s [memory] is/are not supported by SgeBackend. Unsupported attributes will not be part of jobs executions."), occurrences = 1) intercept {
          backend ! Initialize
        }
      }
    }
  }
}
