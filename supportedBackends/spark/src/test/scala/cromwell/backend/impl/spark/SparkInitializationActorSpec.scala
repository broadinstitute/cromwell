package cromwell.backend.impl.spark

import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.TestKitSuite
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s._

import scala.concurrent.duration._

class SparkInitializationActorSpec  extends  TestKitSuite("SparkInitializationActorSpec")
  with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender {

  val Timeout = 5.second.dilated

  val HelloWorld =
    """
      |task hello {
      |  command {
      |    helloApp
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

  private def getSparkBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall], conf: BackendConfigurationDescriptor) = {
    system.actorOf(SparkInitializationActor.props(workflowDescriptor, calls, conf, emptyActor))
  }

  "SparkInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        EventFilter.warning(message = s"Key/s [memory] is/are not supported by SparkBackend. Unsupported attributes will not be part of jobs executions.", occurrences = 1) intercept {
          val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 %s: "%s"}""".format("appMainClass", "test"))
          val backend = getSparkBackend(workflowDescriptor, workflowDescriptor.workflow.taskCalls, emptyBackendConfig)
          backend ! Initialize
        }
      }
    }
  }
}
