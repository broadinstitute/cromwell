package cromwell.backend.impl.spark

import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, TestConfig}
import cromwell.core.TestKitSuite
import org.scalatest.BeforeAndAfterAll
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class SparkInitializationActorSpec  extends  TestKitSuite("SparkInitializationActorSpec")
  with AnyWordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender {

  val Timeout = 10.second.dilated

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

  private def getSparkBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[CommandCallNode], conf: BackendConfigurationDescriptor) = {
    system.actorOf(SparkInitializationActor.props(workflowDescriptor, calls, conf, emptyActor))
  }

  "SparkInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        EventFilter.warning(message = s"Key/s [memory] is/are not supported by SparkBackend. Unsupported attributes will not be part of jobs executions.", occurrences = 1) intercept {
          val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld, runtime = """runtime { memory: 1 %s: "%s"}""".format("appMainClass", "test"))
          val backend = getSparkBackend(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, TestConfig.emptyBackendConfigDescriptor)
          backend ! Initialize
        }
      }
    }
  }
}
