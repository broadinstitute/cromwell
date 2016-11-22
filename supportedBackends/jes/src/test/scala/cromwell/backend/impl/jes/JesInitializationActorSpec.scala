package cromwell.backend.impl.jes

import java.util.UUID

import akka.testkit._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, InitializationSuccess, Initialize}
import cromwell.backend.impl.jes.authentication.GcsLocalizing
import cromwell.backend.{BackendConfigurationDescriptor, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.Tags.IntegrationTest
import cromwell.core.logging.LoggingTest._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.gcs.GoogleConfiguration
import cromwell.filesystems.gcs.auth.{RefreshTokenMode, SimpleClientSecrets}
import cromwell.util.{EncryptionSpec, SampleWdl}
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import spray.json._
import wdl4s.TaskCall

import scala.concurrent.duration._

class JesInitializationActorSpec extends TestKitSuite("JesInitializationActorSpec") with FlatSpecLike with Matchers
  with ImplicitSender with Mockito {
  val Timeout = 5.second.dilated

  import BackendSpec._

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
      |    },
      |    {
      |      name = "user-via-refresh"
      |      scheme = "refresh_token"
      |      client-id = "secret_id"
      |      client-secret = "secret_secret"
      |    }
      |  ]
      |}
      | """.stripMargin)

  val backendConfigTemplate =
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
      |
      |[DOCKERHUBCONFIG]
      |""".stripMargin

  val refreshTokenConfigTemplate =
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
      |      auth = "user-via-refresh"
      |    }
      |  }
      |""".stripMargin

  val backendConfig = ConfigFactory.parseString(backendConfigTemplate.replace("[DOCKERHUBCONFIG]", ""))

  val dockerBackendConfig = ConfigFactory.parseString(backendConfigTemplate.replace("[DOCKERHUBCONFIG]",
    """
      |dockerhub {
      |  account = "my@docker.account"
      |  token = "mydockertoken"
      |}
      | """.stripMargin))

  val defaultBackendConfig = BackendConfigurationDescriptor(backendConfig, globalConfig)

  val refreshTokenConfig = ConfigFactory.parseString(refreshTokenConfigTemplate)

  private def getJesBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall], conf: BackendConfigurationDescriptor) = {
    system.actorOf(JesInitializationActor.props(workflowDescriptor, calls, new JesConfiguration(conf), emptyActor))
  }

  behavior of "JesInitializationActor"

  it should "log a warning message when there are unsupported runtime attributes" taggedAs IntegrationTest in {
    within(Timeout) {
      val workflowDescriptor = buildWorkflowDescriptor(HelloWorld,
        runtime = """runtime { docker: "ubuntu/latest" test: true }""")
      val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflow.taskCalls,
        defaultBackendConfig)
      val eventPattern =
        "Key/s [test] is/are not supported by JesBackend. Unsupported attributes will not be part of jobs executions."
      EventFilter.warning(pattern = escapePattern(eventPattern), occurrences = 1) intercept {
        backend ! Initialize
      }
      expectMsgPF() {
        case InitializationSuccess(_) => //Docker entry is present.
        case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
      }
    }
  }

  it should "return InitializationFailed when docker runtime attribute key is not present" in {
    within(Timeout) {
      val workflowDescriptor = buildWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
      val backend = getJesBackend(workflowDescriptor, workflowDescriptor.workflow.taskCalls,
        defaultBackendConfig)
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

  private case class TestingBits(actorRef: TestActorRef[JesInitializationActor], jesConfiguration: JesConfiguration)

  private def buildJesInitializationTestingBits(backendConfig: Config = dockerBackendConfig): TestingBits = {
    val workflowOptions = WorkflowOptions.fromMap(Map("refresh_token" -> "mytoken")).get
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource(), options = workflowOptions)
    val calls = workflowDescriptor.workflow.taskCalls
    val backendConfigurationDescriptor = BackendConfigurationDescriptor(backendConfig, globalConfig)
    val jesConfiguration = new JesConfiguration(backendConfigurationDescriptor)

    val actorRef = TestActorRef[JesInitializationActor](
      JesInitializationActor.props(workflowDescriptor, calls, jesConfiguration, emptyActor),
      "TestableJesInitializationActor-" + UUID.randomUUID)
    TestingBits(actorRef, jesConfiguration)
  }

  it should "create a GcsLocalizing instance" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, _) = buildJesInitializationTestingBits(refreshTokenConfig)
    val actor = actorRef.underlyingActor
    actor.refreshTokenAuth should be(Some(GcsLocalizing(RefreshTokenMode("user-via-refresh", "secret_id", "secret_secret", GoogleConfiguration.GoogleScopes),  "mytoken")))
  }

  it should "generate the correct json content for no docker token and no refresh token" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, _) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    actor.generateAuthJson(None, None) should be(empty)

    val authJsonOption = actor.generateAuthJson(None, None)
    authJsonOption should be(empty)

    actorRef.stop()
  }

  it should "generate the correct json content for a docker token and no refresh token" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, jesConfiguration) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    val authJsonOption = actor.generateAuthJson(jesConfiguration.dockerCredentials, None)
    authJsonOption shouldNot be(empty)
    authJsonOption.get should be(
      normalize(
        """
          |{
          |    "auths": {
          |        "docker": {
          |            "account": "my@docker.account",
          |            "token": "mydockertoken"
          |        }
          |    }
          |}
        """.stripMargin)
    )

    actorRef.stop()
  }

  it should "generate the correct json content for no docker token and a refresh token" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, _) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    val gcsUserAuth = Option(GcsLocalizing(SimpleClientSecrets("myclientid", "myclientsecret"), "mytoken"))
    val authJsonOption = actor.generateAuthJson(None, gcsUserAuth)
    authJsonOption shouldNot be(empty)
    authJsonOption.get should be(
      normalize(
        """
          |{
          |    "auths": {
          |        "boto": {
          |            "client_id": "myclientid",
          |            "client_secret": "myclientsecret",
          |            "refresh_token": "mytoken"
          |        }
          |    }
          |}
        """.stripMargin)
    )

    actorRef.stop()
  }

  it should "generate the correct json content for a docker token and a refresh token" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, jesConfiguration) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    val gcsUserAuth = Option(GcsLocalizing(SimpleClientSecrets("myclientid", "myclientsecret"), "mytoken"))
    val authJsonOption = actor.generateAuthJson(jesConfiguration.dockerCredentials, gcsUserAuth)
    authJsonOption shouldNot be(empty)
    authJsonOption.get should be(
      normalize(
        """
          |{
          |    "auths": {
          |        "docker": {
          |            "account": "my@docker.account",
          |            "token": "mydockertoken"
          |        },
          |        "boto": {
          |            "client_id": "myclientid",
          |            "client_secret": "myclientsecret",
          |            "refresh_token": "mytoken"
          |        }
          |    }
          |}
        """.stripMargin)
    )

    actorRef.stop()
  }

  private def normalize(str: String) = {
    str.parseJson.prettyPrint
  }
}
