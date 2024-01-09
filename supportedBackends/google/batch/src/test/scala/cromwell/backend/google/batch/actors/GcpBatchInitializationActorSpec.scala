package cromwell.backend.google.batch.actors

import java.util.UUID
import akka.actor.Props
import akka.testkit._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.async.RuntimeAttributeValidationFailures
import cromwell.backend.google.batch.models.GcpBatchConfiguration
import cromwell.backend.google.batch.actors.GcpBatchInitializationActorSpec._
import cromwell.backend.google.batch.models.GcpBatchTestConfig.{batchAttributes, googleConfiguration, BatchGlobalConfig}
import cromwell.backend.{BackendConfigurationDescriptor, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.TestKitSuite
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.logging.LoggingTest._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class GcpBatchInitializationActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with ImplicitSender {
  val Timeout: FiniteDuration = 30.second.dilated

  import BackendSpec._

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

  private def getJesBackendProps(workflowDescriptor: BackendWorkflowDescriptor,
                                 calls: Set[CommandCallNode],
                                 jesConfiguration: GcpBatchConfiguration
  ): Props = {
    val ioActor = mockIoActor
    val params = GcpBatchInitializationActorParams(workflowDescriptor,
                                                   ioActor,
                                                   calls,
                                                   jesConfiguration,
                                                   emptyActor,
                                                   restarting = false
    )
    Props(new GcpBatchInitializationActor(params)).withDispatcher(BackendDispatcher)
  }

  private def getJesBackend(workflowDescriptor: BackendWorkflowDescriptor,
                            calls: Set[CommandCallNode],
                            conf: BackendConfigurationDescriptor
  ) = {
    val props = getJesBackendProps(workflowDescriptor,
                                   calls,
                                   new GcpBatchConfiguration(conf, googleConfiguration, batchAttributes)
    )
    system.actorOf(props, "TestableJesInitializationActor-" + UUID.randomUUID)
  }

  behavior of "GcpBatchInitializationActor"

  it should "log a warning message when there are unsupported runtime attributes" in {

    within(Timeout) {
      val workflowDescriptor =
        buildWdlWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu/latest" test: true }""")
      val backend = getJesBackend(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, defaultBackendConfig)
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

  it should "return InitializationFailed when docker runtime attribute key is not present" in {
    within(Timeout) {
      val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
      val backend = getJesBackend(workflowDescriptor, workflowDescriptor.callable.taskCallNodes, defaultBackendConfig)
      backend ! Initialize
      expectMsgPF() { case InitializationFailed(failure) =>
        failure match {
          case exception: RuntimeAttributeValidationFailures =>
            if (
              !exception.getMessage.equals(
                "Runtime validation failed:\nTask hello has an invalid runtime attribute docker = !! NOT FOUND !!"
              )
            )
              fail(
                "Exception message is not equal to 'Runtime validation failed:\nTask hello has an invalid runtime attribute docker = !! NOT FOUND !!'."
              )
        }
      }
    }
  }
}

object GcpBatchInitializationActorSpec {
  val globalConfig: Config = ConfigFactory.parseString("""
                                                         |google {
                                                         |
                                                         |  application-name = "cromwell"
                                                         |
                                                         |  auths = [
                                                         |    {
                                                         |      name = "application-default"
                                                         |      scheme = "mock"
                                                         |    }
                                                         |  ]
                                                         |}
                                                         |""".stripMargin)

  val backendConfigTemplate: String =
    """
      |  // Google project
      |  project = "my-cromwell-workflows"
      |
      |  // Base bucket for workflow executions
      |  root = "gs://my-cromwell-workflows-bucket"
      |
      |  // Polling for completion backs-off gradually for slower-running jobs.
      |  // This is the maximum polling interval (in seconds):
      |  maximum-polling-interval = 600
      |
      |  genomics {
      |  // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
      |  // Pipelines and manipulate auth JSONs.
      |     auth = "application-default"
      |     // Endpoint for APIs, no reason to change this unless directed by Google.
      |     endpoint-url = "https://genomics.googleapis.com/"
      |  }
      |
      |  default-runtime-attributes {
      |     cpu: 1
      |     failOnStderr: false
      |     # Allowed to be a boolean, or a list of Ints, or an Int
      |     continueOnReturnCode: 0
      |     memory: "2 GB"
      |     bootDiskSizeGb: 10
      |     # Allowed to be a String, or a list of Strings
      |     disks: "local-disk 10 SSD"
      |     noAddress: false
      |     preemptible: 0
      |     zones: ["us-central1-a", "us-central1-b"]
      |  }
      |  filesystems {
      |    gcs {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "application-default"
      |    }
      |  }
      |
      |[VPCCONFIG]
      |
      |[DOCKERHUBCONFIG]
      |""".stripMargin

  val backendConfig: Config =
    ConfigFactory.parseString(backendConfigTemplate.replace("[VPCCONFIG]", "").replace("[DOCKERHUBCONFIG]", ""))

  val dockerBackendConfig: Config = ConfigFactory.parseString(
    backendConfigTemplate
      .replace("[VPCCONFIG]", "")
      .replace(
        "[DOCKERHUBCONFIG]",
        """
          |dockerhub {
          |  account = "my@docker.account"
          |  # no secrets here guys this is just `echo -n username:password | base64`
          |  token = "dXNlcm5hbWU6cGFzc3dvcmQ="
          |}
          | """.stripMargin
      )
  )

  val vpcBackendConfig: Config = ConfigFactory.parseString(
    backendConfigTemplate
      .replace("[DOCKERHUBCONFIG]", "")
      .replace(
        "[VPCCONFIG]",
        """
          |virtual-private-cloud {
          |  network-label-key = "cromwell-ci-network"
          |  subnetwork-label-key = "cromwell-ci-subnetwork"
          |  auth = "service_account"
          |}
          | """.stripMargin
      )
  )

  private val defaultBackendConfig = new BackendConfigurationDescriptor(backendConfig, globalConfig) {
    override private[backend] lazy val cromwellFileSystems = new CromwellFileSystems(BatchGlobalConfig)
  }
}
