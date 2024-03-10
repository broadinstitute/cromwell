package cromwell.backend.google.batch.models

import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.batch.models.GcpBatchTestConfig._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class GcpBatchConfigurationSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks
    with BeforeAndAfterAll {

  behavior of "GcpBatchConfigurationSpec"

  val mockFile = DefaultPathBuilder.createTempFile()

  override def afterAll(): Unit = {
    mockFile.delete(swallowIOExceptions = true)
    ()
  }

  val globalConfig = ConfigFactory.parseString(s"""
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
                                                  |      name = "service-account"
                                                  |      scheme = "service_account"
                                                  |      service-account-id = "my-service-account"
                                                  |      pem-file = "${mockFile.pathAsString}"
                                                  |    }
                                                  |  ]
                                                  |}
                                                  |
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
      |     location = "us-central1"
      |  }
      |
      |  default-runtime-attributes {
      |      failOnStderr: false
      |      continueOnReturnCode: 0
      |      cpu: 1
      |      memory: "2 GB"
      |      bootDiskSizeGb: 10
      |      disks: "local-disk 10 SSD"
      |      noAddress: false
      |      preemptible: 3
      |      zones:["us-central1-a", "us-central1-b"]
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
      |
    """.stripMargin
  )

  it should "fail to instantiate if any required configuration is missing" in {

    val configs = Table(
      ("backendConfig", "globalConfig"),
      (backendConfig, globalConfig.withoutPath("google")),
      (backendConfig.withoutPath("project"), globalConfig),
      (backendConfig.withoutPath("root"), globalConfig),
      (backendConfig.withoutPath("genomics"), globalConfig),
      (backendConfig.withoutPath("genomics.location"), globalConfig),
      (backendConfig.withoutPath("filesystems"), globalConfig),
      (backendConfig.withoutPath("filesystems.gcs"), globalConfig),
      (backendConfig.withoutPath("filesystems.gcs.auth"), globalConfig)
    )

    forAll(configs) { (backend, global) =>
      an[Exception] shouldBe thrownBy {
        val failingGoogleConf = GoogleConfiguration(global)
        val failingAttributes = GcpBatchConfigurationAttributes(failingGoogleConf, backend, "papi")
        new GcpBatchConfiguration(BackendConfigurationDescriptor(backend, global), failingGoogleConf, failingAttributes)
      }
    }
  }

  it should "have correct root" in {
    new GcpBatchConfiguration(BackendConfigurationDescriptor(backendConfig, globalConfig),
                              googleConfiguration,
                              batchAttributes
    ).root shouldBe "gs://my-cromwell-workflows-bucket"
  }

  it should "have correct docker" in {
    val dockerConf = new GcpBatchConfiguration(
      BackendConfigurationDescriptor(backendConfig, globalConfig),
      googleConfiguration,
      batchAttributes
    ).dockerCredentials

    dockerConf shouldBe defined
    dockerConf.get.token shouldBe "dockerToken"
  }
}
