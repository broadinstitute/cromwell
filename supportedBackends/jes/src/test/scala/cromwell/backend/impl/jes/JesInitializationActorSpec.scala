package cromwell.backend.impl.jes

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import cromwell.core.logging.LoggingTest._
import spray.json.{JsObject, JsValue}
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlSource}

import scala.concurrent.duration._

class JesInitializationActorSpec extends TestKit(ActorSystem("JesInitializationActorSpec", ConfigFactory.parseString(
  // TODO: PBE: 5s leeway copy of CromwellTestkitSpec. Refactor to D.R.Y. this code, and rename Testkit to TestKit
  """
    |akka {
    |  loggers = ["akka.testkit.TestEventListener"]
    |  loglevel = "INFO"
    |  actor {
    |    debug {
    |       receive = on
    |    }
    |  }
    |  dispatchers {
    |    slow-actor-dispatcher {
    |      type = Dispatcher
    |      executor = "fork-join-executor"
    |    }
    |  }
    |  test {
    |    # Some of our tests fire off a message, then expect a particular event message within 3s (the default).
    |    # Especially on CI, the metadata test does not seem to be returning in time. So, overriding the timeouts
    |    # with slightly higher values. Alternatively, could also adjust the akka.test.timefactor only in CI.
    |    filter-leeway = 5s
    |    single-expect-default = 5s
    |    default-timeout = 10s
    |  }
    |}
    |""".stripMargin))) with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender {
  val Timeout = 5.second.dilated

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

  val globalConfig = ConfigFactory.parseString(
    """
      |google {
      |
      |  application-name = "cromwell"
      |
      |  auths = [
      |    {
      |      name = "application-default"
      |      scheme = "application_default"
      |    }
      |  ]
      |}
    """.stripMargin)
  val backendConfig = ConfigFactory.parseString(
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
      |  filesystems {
      |    gcs {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "application-default"
      |    }
      |  }
    """.stripMargin)

  val defaultBackendConfig = new BackendConfigurationDescriptor(backendConfig, globalConfig)

  private def buildWorkflowDescriptor(wdl: WdlSource,
                                      inputs: Map[String, WdlValue] = Map.empty,
                                      options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                                      runtime: String = "") = {
    new BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime)),
      inputs,
      options
    )
  }

  private def getJesBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], conf: BackendConfigurationDescriptor) = {
    system.actorOf(JesInitializationActor.props(workflowDescriptor, calls, new JesConfiguration(conf)))
  }

  override def afterAll {
    system.shutdown()
  }

  "JesInitializationActor" should {
    "log a warning message when there are unsupported runtime attributes" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { docker: "ubuntu/latest" test: true }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        EventFilter.warning(pattern = escapePattern(s"Key/s [test] is/are not supported by JesBackend. Unsupported attributes will not be part of jobs executions."), occurrences = 1) intercept {
          backend ! Initialize
        }
        expectMsgPF() {
          case InitializationSuccess => //Docker entry is present.
          case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
        }
      }
    }

    "return InitializationFailed when docker runtime attribute key is not present" in {
      within(Timeout) {
        val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
        val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflowNamespace.workflow.calls, defaultBackendConfig)
        backend ! Initialize
        expectMsgPF() {
          case InitializationFailed(failure) =>
            failure match {
              case exception: IllegalArgumentException =>
                if (!exception.getMessage.equals("Task hello has an invalid runtime attribute docker = !! NOT FOUND !!"))
                  fail("Exception message does not contains 'Runtime attribute validation failed'.")
            }
        }
      }
    }
  }
}
