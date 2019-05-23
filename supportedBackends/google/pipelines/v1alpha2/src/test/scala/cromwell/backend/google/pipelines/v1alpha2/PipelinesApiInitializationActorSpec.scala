package cromwell.backend.google.pipelines.v1alpha2

import java.util.UUID

import akka.actor.Props
import akka.testkit._
import com.typesafe.config.{Config, ConfigFactory}
import common.collections.EnhancedCollections._
import cromwell.backend.google.pipelines.common.PipelinesApiInitializationActorSpec._
import cromwell.backend.google.pipelines.common.PipelinesApiTestConfig.genomicsFactory
import cromwell.backend.google.pipelines.common.authentication.{GcsLocalizing, PipelinesApiAuthObject}
import cromwell.backend.google.pipelines.common.{PipelinesApiConfigurationAttributes, PipelinesApiConfiguration, PipelinesApiInitializationActorParams}
import cromwell.backend.{BackendConfigurationDescriptor, BackendSpec, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.{RefreshTokenMode, SimpleClientSecrets}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.util.{EncryptionSpec, SampleWdl}
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import spray.json._
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class PipelinesApiInitializationActorSpec extends TestKitSuite("PipelinesApiInitializationActorSpecV1") with FlatSpecLike with Matchers
  with ImplicitSender with Mockito {
  val Timeout: FiniteDuration = 10.second.dilated

  import BackendSpec._
  import PipelinesApiInitializationActorSpec._

  val refreshTokenConfigTemplate: String =
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
      |
      |  filesystems {
      |    gcs {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "user-via-refresh"
      |    }
      |  }
      |""".stripMargin

  val refreshTokenConfig: Config = ConfigFactory.parseString(refreshTokenConfigTemplate)

  private def getJesBackendProps(workflowDescriptor: BackendWorkflowDescriptor,
                                 calls: Set[CommandCallNode],
                                 jesConfiguration: PipelinesApiConfiguration): Props = {
    val ioActor = mockIoActor
    val params = PipelinesApiInitializationActorParams(workflowDescriptor, ioActor, calls, jesConfiguration, emptyActor, restarting = false)
    Props(new PipelinesApiInitializationActor(params)).withDispatcher(BackendDispatcher)
  }

  behavior of "PipelinesApiInitializationActorV1"

  private case class TestingBits(actorRef: TestActorRef[PipelinesApiInitializationActor], jesConfiguration: PipelinesApiConfiguration)

  private def buildJesInitializationTestingBits(backendConfig: Config = dockerBackendConfig): TestingBits = {
    val workflowOptions = WorkflowOptions.fromMap(Map("refresh_token" -> "mytoken")).get
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint),
      options = workflowOptions
    )
    val calls = workflowDescriptor.callable.taskCallNodes
    val backendConfigurationDescriptor = BackendConfigurationDescriptor(backendConfig, globalConfig)
    val customGoogleConfig = GoogleConfiguration(globalConfig)
    val customAttributes = PipelinesApiConfigurationAttributes(customGoogleConfig, backendConfig)
    val jesConfiguration = new PipelinesApiConfiguration(backendConfigurationDescriptor, genomicsFactory, customGoogleConfig, customAttributes)

    val actorRef = TestActorRef[PipelinesApiInitializationActor](
      getJesBackendProps(workflowDescriptor, calls, jesConfiguration),
      "TestableJesInitializationActor-" + UUID.randomUUID)
    TestingBits(actorRef, jesConfiguration)
  }

  it should "create a GcsLocalizing instance" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, _) = buildJesInitializationTestingBits(refreshTokenConfig)
    val actor = actorRef.underlyingActor
    val expectedAuthMode = RefreshTokenMode("user-via-refresh", "secret_id", "secret_secret")
    val expectedAuth = GcsLocalizing(expectedAuthMode,  "mytoken")
    actor.refreshTokenAuth should be(Option(expectedAuth))
  }

  it should "generate the correct json content for no docker token and no refresh token" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, _) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    actor.generateAuthJson(flattenAuthOptions(None, None), restrictMetadataAccess = false) should be(empty)

    val authJsonOption = actor.generateAuthJson(flattenAuthOptions(None, None), restrictMetadataAccess = false)
    authJsonOption should be(empty)

    actorRef.stop()
  }

  it should "generate the correct json content for a docker token and no refresh token" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, jesConfiguration) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    val authJsonOption = actor.generateAuthJson(
      flattenAuthOptions(jesConfiguration.dockerCredentials, None),
      restrictMetadataAccess = false
    )
    authJsonOption shouldNot be(empty)
    // dXNlcm5hbWU6cGFzc3dvcmQ= is base64-encoded username:password
    authJsonOption.get should be(
      normalize(
        """
          |{
          |    "auths": {
          |        "docker": {
          |            "token": "dXNlcm5hbWU6cGFzc3dvcmQ="
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
    val authJsonOption = actor.generateAuthJson(flattenAuthOptions(None, gcsUserAuth), restrictMetadataAccess = false)
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
    val authJsonOption = actor.generateAuthJson(
      flattenAuthOptions(jesConfiguration.dockerCredentials, gcsUserAuth),
      restrictMetadataAccess = false
    )
    authJsonOption shouldNot be(empty)
    normalize(authJsonOption.get) should be(
      normalize(
        """
          |{
          |    "auths": {
          |        "docker": {
          |            "token": "dXNlcm5hbWU6cGFzc3dvcmQ="
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

  it should "generate the correct json content for a docker token, a refresh token, and restrictMetadataAccess" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, jesConfiguration) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    val gcsUserAuth = Option(GcsLocalizing(SimpleClientSecrets("myclientid", "myclientsecret"), "mytoken"))
    val authJsonOption = actor.generateAuthJson(
      flattenAuthOptions(jesConfiguration.dockerCredentials, gcsUserAuth),
      restrictMetadataAccess = true
    )
    authJsonOption shouldNot be(empty)
    normalize(authJsonOption.get) should be(
      normalize(
        """
          |{
          |    "auths": {
          |        "docker": {
          |            "token": "dXNlcm5hbWU6cGFzc3dvcmQ="
          |        },
          |        "boto": {
          |            "client_id": "myclientid",
          |            "client_secret": "myclientsecret",
          |            "refresh_token": "mytoken"
          |        }
          |    },
          |    "restrictMetadataAccess": true
          |}
        """.stripMargin)
    )

    actorRef.stop()
  }

  it should "generate the correct json content for just restrictMetadataAccess" in {
    EncryptionSpec.assumeAes256Cbc()

    val TestingBits(actorRef, _) = buildJesInitializationTestingBits()
    val actor = actorRef.underlyingActor

    val authJsonOption = actor.generateAuthJson(flattenAuthOptions(None, None), restrictMetadataAccess = true)
    authJsonOption shouldNot be(empty)
    authJsonOption.get should be(
      normalize(
        """
          |{
          |    "restrictMetadataAccess": true
          |}
        """.stripMargin)
    )

    actorRef.stop()
  }
}

object PipelinesApiInitializationActorSpec {
  def normalize(str: String) = {
    str.parseJson.prettyPrint
  }

  def flattenAuthOptions(options: Option[PipelinesApiAuthObject]*): List[PipelinesApiAuthObject] = {
    options.toList.flatten
  }
}
