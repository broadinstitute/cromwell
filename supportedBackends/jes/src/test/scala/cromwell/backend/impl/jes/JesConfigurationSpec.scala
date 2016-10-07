package cromwell.backend.impl.jes

import better.files.File
import com.typesafe.config.{ConfigValueFactory, ConfigFactory}
import cromwell.backend.BackendConfigurationDescriptor
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

class JesConfigurationSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with BeforeAndAfterAll {

  behavior of "JesConfigurationSpec"

  val mockFile = File.newTemporaryFile()

  override def afterAll(): Unit = {
    mockFile.delete(true)

    ()
  }

  val globalConfig = ConfigFactory.parseString(
    s"""
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
      |      client-secret = "${mockFile.pathAsString}"
      |    },
      |    {
      |      name = "service-account"
      |      scheme = "service_account"
      |      service-account-id = "my-service-account"
      |      pem-file = "${mockFile.pathAsString}"
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
      |  dockerhub {
      |    account = "dockerAccount"
      |    token = "dockerToken"
      |  }
      |
      |  filesystems {
      |    gcs {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "application-default"
      |    }
      |  }
    """.stripMargin)

  it should "fail to instantiate if any required configuration is missing" in {

    val configs = Table(
      ("backendConfig", "globalConfig"),
      (backendConfig, globalConfig.withoutPath("google")),
      (backendConfig.withoutPath("project"), globalConfig),
      (backendConfig.withoutPath("root"), globalConfig),
      (backendConfig.withoutPath("genomics"), globalConfig),
      (backendConfig.withoutPath("genomics.endpoint-url"), globalConfig),
      (backendConfig.withoutPath("filesystems"), globalConfig),
      (backendConfig.withoutPath("filesystems.gcs"), globalConfig),
      (backendConfig.withoutPath("filesystems.gcs.auth"), globalConfig)
    )

    forAll(configs) { (backend, global) =>
      an[Exception] shouldBe thrownBy {
        new JesConfiguration(BackendConfigurationDescriptor(backend, global))
      }
    }
  }

  it should "have correct root" in {
    new JesConfiguration(BackendConfigurationDescriptor(backendConfig, globalConfig)).root shouldBe "gs://my-cromwell-workflows-bucket"
  }

  it should "have correct docker" in {
    val dockerConf = new JesConfiguration(BackendConfigurationDescriptor(backendConfig, globalConfig)).dockerCredentials
    dockerConf shouldBe defined
    dockerConf.get.account shouldBe "dockerAccount"
    dockerConf.get.token shouldBe "dockerToken"
  }

  it should "have correct needAuthFileUpload" in {
    val configs = Table(
      ("backendConfig", "globalConfig"),
      // With Docker
      (backendConfig.withValue("filesystems.gcs.auth", ConfigValueFactory.fromAnyRef("application-default")), true),
      (backendConfig.withValue("filesystems.gcs.auth", ConfigValueFactory.fromAnyRef("user-via-refresh")), true),
      (backendConfig.withValue("filesystems.gcs.auth", ConfigValueFactory.fromAnyRef("service-account")), true),
      // Without Docker
      (backendConfig.withValue("filesystems.gcs.auth", ConfigValueFactory.fromAnyRef("application-default")).withoutPath("dockerhub"), false),
      (backendConfig.withValue("filesystems.gcs.auth", ConfigValueFactory.fromAnyRef("user-via-refresh")).withoutPath("dockerhub"), true),
      (backendConfig.withValue("filesystems.gcs.auth", ConfigValueFactory.fromAnyRef("service-account")).withoutPath("dockerhub"), false)
    )

    forAll(configs) { (backend, needAuthFileUpload) =>
      new JesConfiguration(BackendConfigurationDescriptor(backend, globalConfig)).needAuthFileUpload shouldBe needAuthFileUpload
    }
  }
}
