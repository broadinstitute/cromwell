package cromwell.backend.sfs

import akka.actor.{ActorRef, Props}
import akka.testkit.{EventFilter, ImplicitSender, TestDuration}
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import cromwell.backend.BackendSpec._
import cromwell.backend.BackendWorkflowInitializationActor.Initialize
import cromwell.backend.standard.DefaultInitializationActorParams
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, TestConfig}
import cromwell.core.TestKitSuite
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.logging.LoggingTest._
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class SharedFileSystemInitializationActorSpec extends TestKitSuite("SharedFileSystemInitializationActorSpec")
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

  private def getActorRef(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[CommandCallNode],
                          conf: BackendConfigurationDescriptor) = {
    val params = DefaultInitializationActorParams(workflowDescriptor, emptyActor, calls, emptyActor, conf, restarting = false)
    val props = Props(new SharedFileSystemInitializationActor(params))
    system.actorOf(props, "SharedFileSystemInitializationActor")
  }

  "SharedFileSystemInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld, runtime = """runtime { unsupported: 1 }""")
        val mockFileSystems = new CromwellFileSystems(ConfigFactory.empty())
        val conf = new BackendConfigurationDescriptor(TestConfig.sampleBackendRuntimeConfig, ConfigFactory.empty()) {
          override lazy val configuredPathBuilderFactories = mockFileSystems.factoriesFromConfig(TestConfig.sampleBackendRuntimeConfig).unsafe("Failed to instantiate backend filesystem")
        }
        val backend: ActorRef = getActorRef(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, conf)
        val pattern = "Key/s [unsupported] is/are not supported by backend. " +
          "Unsupported attributes will not be part of job executions."
        EventFilter.warning(pattern = escapePattern(pattern), occurrences = 1) intercept {
          backend ! Initialize
        }
      }
    }
  }
}
