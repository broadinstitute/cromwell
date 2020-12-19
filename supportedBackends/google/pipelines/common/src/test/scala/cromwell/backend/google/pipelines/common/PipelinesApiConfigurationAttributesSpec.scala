package cromwell.backend.google.pipelines.common

import java.net.URL

import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.{Config, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import common.exception.MessageAggregation
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.filesystems.gcs.GcsPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._

class PipelinesApiConfigurationAttributesSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  import PipelinesApiTestConfig._

  behavior of "PipelinesApiAttributes"

  val googleConfig = GoogleConfiguration(PapiGlobalConfig)
  val runtimeConfig: Config = ConfigFactory.load()

  it should "parse correct PAPI config" in {

    val backendConfig = ConfigFactory.parseString(configString())

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.endpointUrl should be(new URL("http://myEndpoint"))
    pipelinesApiAttributes.project should be("myProject")
    pipelinesApiAttributes.executionBucket should be("gs://myBucket")
    pipelinesApiAttributes.maxPollingInterval should be(600)
    pipelinesApiAttributes.computeServiceAccount should be("default")
    pipelinesApiAttributes.restrictMetadataAccess should be(false)
    pipelinesApiAttributes.memoryRetryConfiguration should be(None)
    pipelinesApiAttributes.referenceFileToDiskImageMappingOpt.isEmpty should be(true)
  }

  it should "parse correct preemptible config" in {
    val backendConfig = ConfigFactory.parseString(configString(customContent = "preemptible = 3"))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.endpointUrl should be(new URL("http://myEndpoint"))
    pipelinesApiAttributes.project should be("myProject")
    pipelinesApiAttributes.executionBucket should be("gs://myBucket")
    pipelinesApiAttributes.maxPollingInterval should be(600)
  }

  it should "parse batch-requests.timeouts values correctly" in {
    val customContent =
      """
        |batch-requests {
        |  timeouts {
        |    read = 100 hours
        |    connect = 10 seconds
        |  }
        |}
      """.stripMargin

    val backendConfig = ConfigFactory.parseString(configString(customContent = customContent))
    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")

    pipelinesApiAttributes.batchRequestTimeoutConfiguration.readTimeoutMillis.get.value should be(100.hours.toMillis.toInt)
    pipelinesApiAttributes.batchRequestTimeoutConfiguration.connectTimeoutMillis.get.value should be(10.seconds.toMillis.toInt)
  }

  it should "parse an empty batch-requests.timeouts section correctly" in {
    val customContent =
      """
        |batch-requests {
        |  timeouts {
        |    # read = 100 hours
        |    # connect = 10 seconds
        |  }
        |}
      """.stripMargin

    val backendConfig = ConfigFactory.parseString(configString(customContent = customContent))
    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")

    pipelinesApiAttributes.batchRequestTimeoutConfiguration should be(BatchRequestTimeoutConfiguration(None, None))
  }

  it should "parse pipeline-timeout" in {
    val backendConfig = ConfigFactory.parseString(configString(customContent = "pipeline-timeout = 3 days"))
    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")

    pipelinesApiAttributes.pipelineTimeout should be(3.days)
  }

  it should "parse an undefined pipeline-timeout" in {
    val backendConfig = ConfigFactory.parseString(configString())
    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")

    pipelinesApiAttributes.pipelineTimeout should be(7.days)
  }

  it should "parse compute service account" in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = """compute-service-account = "testing" """))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.computeServiceAccount should be("testing")
  }

  it should "parse restrict-metadata-access" in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = "restrict-metadata-access = true"))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.restrictMetadataAccess should be(true)

  }

  it should "parse localization-attempts" in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = "localization-attempts = 31380"))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.gcsTransferConfiguration.transferAttempts.value should be(31380)
  }

  it should "parse virtual-private-cloud" in {

    val customConfig =
      """
        |  virtual-private-cloud {
        |    network-label-key = my-network
        |    subnetwork-label-key = my-subnetwork
        |    auth = application-default
        |  }
      """.stripMargin

    val backendConfig = ConfigFactory.parseString(configString(customConfig))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.name should be("my-network")
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.subnetwork should be (Option("my-subnetwork"))
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.auth.name should be("application-default")
  }

  it should "parse virtual-private-cloud without subnetwork key" in {

    val customConfig =
      """
        |  virtual-private-cloud {
        |    network-label-key = my-network
        |    auth = application-default
        |  }
      """.stripMargin

    val backendConfig = ConfigFactory.parseString(configString(customConfig))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.name should be("my-network")
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.subnetwork should be(None)
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.auth.name should be("application-default")
  }

  it should "parse virtual-private-cloud with empty body" in {

    val customConfig = """virtual-private-cloud{ }"""

    val backendConfig = ConfigFactory.parseString(configString(customConfig))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.virtualPrivateCloudConfiguration should be(None)
  }

  it should "parse memory-retry" in {
    val customConfig =
      """
        |memory-retry {
        |   error-keys = ["OutOfMemory", "Killed", "Exit123"]
        |   multiplier = 1.1
        |}""".stripMargin
    val backendConfig = ConfigFactory.parseString(configString(customConfig))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.memoryRetryConfiguration.get.errorKeys shouldBe List("OutOfMemory", "Killed", "Exit123")
    pipelinesApiAttributes.memoryRetryConfiguration.get.multiplier.value shouldBe 1.1
  }

  it should "parse memory-retry with only error-keys" in {
    val customConfig =
      """
        |memory-retry {
        |   error-keys = ["OutOfMemory", "Killed", "Exit123"]
        |}""".stripMargin
    val backendConfig = ConfigFactory.parseString(configString(customConfig))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.memoryRetryConfiguration.get.errorKeys shouldBe List("OutOfMemory", "Killed", "Exit123")
    pipelinesApiAttributes.memoryRetryConfiguration.get.multiplier.value shouldBe 2.0
  }

  it should "parse memory-retry with empty body" in {
    val customConfig = """memory-retry { }""".stripMargin
    val backendConfig = ConfigFactory.parseString(configString(customConfig))

    val pipelinesApiAttributes = PipelinesApiConfigurationAttributes(googleConfig, backendConfig, "papi")
    pipelinesApiAttributes.memoryRetryConfiguration shouldBe None
  }

  it should "not parse invalid config" in {
    val nakedConfig =
      ConfigFactory.parseString(
        """
          |{
          |   genomics {
          |     endpoint-url = "myEndpoint"
          |   }
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, nakedConfig, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("String: 2: No configuration setting found for key 'project'")
    errorsList should contain("String: 2: No configuration setting found for key 'root'")
    errorsList should contain("String: 3: No configuration setting found for key 'genomics.auth'")
    errorsList should contain("String: 2: No configuration setting found for key 'filesystems'")
    errorsList should contain("String: 2: genomics.endpoint-url has type String rather than java.net.URL")
  }

  it should "not parse invalid virtual-private-cloud config without auth" in {
    val networkKeyOnlyConfig =
      ConfigFactory.parseString(
        """
          |{
          |   virtual-private-cloud {
          |     network-label-key = "my-network"
          |   }
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, networkKeyOnlyConfig, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Virtual Private Cloud configuration is invalid. Missing keys: `auth`.")
  }

  it should "not parse invalid virtual-private-cloud config without network label key" in {
    val authOnlyConfig =
      ConfigFactory.parseString(
        """
          |{
          |   virtual-private-cloud {
          |     auth = "application-default"
          |   }
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, authOnlyConfig, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Virtual Private Cloud configuration is invalid. Missing keys: `network-label-key`.")
  }

  it should "not parse invalid virtual-private-cloud config with just subnetwork label key" in {
    val subnetworkOnlyConfig =
      ConfigFactory.parseString(
        """
          |{
          |   virtual-private-cloud {
          |     subnetwork-label-key = "my-subnetwork"
          |   }
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, subnetworkOnlyConfig, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Virtual Private Cloud configuration is invalid. Missing keys: `network-label-key,auth`.")
  }

  it should "not parse invalid virtual-private-cloud config with network & subnetwork label keys" in {
    val config =
      ConfigFactory.parseString(
        """
          |{
          |   virtual-private-cloud {
          |     network-label-key = "my-network"
          |     subnetwork-label-key = "my-subnetwork"
          |   }
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, config, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Virtual Private Cloud configuration is invalid. Missing keys: `auth`.")
  }

  it should "not parse invalid virtual-private-cloud config with subnetwork label key & auth" in {
    val config =
      ConfigFactory.parseString(
        """
          |{
          |   virtual-private-cloud {
          |     subnetwork-label-key = "my-subnetwork"
          |     auth = "application-default"
          |   }
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, config, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Virtual Private Cloud configuration is invalid. Missing keys: `network-label-key`.")
  }

  it should "not parse memory-retry without error-keys" in {
    val config =
      ConfigFactory.parseString(
        """
          |memory-retry {
          |   multiplier = 1.1
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, config, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("memory-retry configuration is invalid. No error-keys provided.")
  }

  it should "not allow a positive multiplier less than 1.0" in {
    val config =
      ConfigFactory.parseString(
        """
          |memory-retry {
          |   error-keys = ["OutOfMemory", "Killed", "Exit123"]
          |   multiplier = 0.5
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, config, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Value 0.5 for memory-retry.multiplier should be greater than 1.0.")
  }

  it should "not allow multiplier negative multiplier" in {
    val config =
      ConfigFactory.parseString(
        """
          |memory-retry {
          |   error-keys = ["OutOfMemory", "Killed", "Exit123"]
          |   multiplier = -2.0
          |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      PipelinesApiConfigurationAttributes(googleConfig, config, "papi")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Value -2.0 for memory-retry.multiplier should be greater than 1.0.")
  }

  def configString(customContent: String = "", genomics: String = ""): String =
    s"""
      |{
      |   project = "myProject"
      |   root = "gs://myBucket"
      |   maximum-polling-interval = 600
      |   $customContent
      |   genomics {
      |     // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
      |     // Pipelines and manipulate auth JSONs.
      |     auth = "application-default"
      |    $genomics
      |     endpoint-url = "http://myEndpoint"
      |   }
      |
      |   filesystems = {
      |     gcs {
      |       // A reference to a potentially different auth for manipulating files via engine functions.
      |       auth = "application-default"
      |     }
      |   }
      |}
      | """.stripMargin

  it should "parse gsutil memory specifications" in {
    val valids = List("0", "150M", "14   PIBIT", "6kib")

    valids foreach {
      case PipelinesApiConfigurationAttributes.GsutilHumanBytes(_, _) =>
      case bad => fail(s"'$bad' was expected to be a valid gsutil memory specification")
    }
  }

  it should "reject invalid memory specifications" in {
    val invalids = List("-1", "150MB", "14PB")

    invalids foreach {
      case invalid@PipelinesApiConfigurationAttributes.GsutilHumanBytes(_, _) => fail(s"Memory specification $invalid not expected to be accepted")
      case _ =>
    }
  }

  it should "parse correct existing reference-disk-localization-manifest-files config" in {
    val referenceDiskManifest1Path = "gs://bucket/manifest1.json"
    val referenceDiskManifest2Path = "gs://bucket/manifest2.json"
    val referenceDiskManifestConfigStr = s"""reference-disk-localization-manifest-files = ["$referenceDiskManifest1Path", "$referenceDiskManifest2Path"]"""
    val backendConfig = ConfigFactory.parseString(configString(referenceDiskManifestConfigStr))

    val validatedGcsPathsToReferenceDiskManifestFilesErrorOr = PipelinesApiConfigurationAttributes.validateGcsPathsToReferenceDiskManifestFiles(backendConfig)
    validatedGcsPathsToReferenceDiskManifestFilesErrorOr match {
      case Valid(validatedGcsPathsToReferenceDiskManifestFilesOpt) =>
        validatedGcsPathsToReferenceDiskManifestFilesOpt match {
          case Some(validatedGcsPathsToReferenceDiskManifestFiles) =>
            validatedGcsPathsToReferenceDiskManifestFiles should contain allElementsOf
              List(
                GcsPathBuilder.validateGcsPath(referenceDiskManifest1Path),
                GcsPathBuilder.validateGcsPath(referenceDiskManifest2Path)
              )
          case None =>
            fail("GCS paths to reference disk manifest files, parsed from config, should not be empty")
        }
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to reference disk manifest files from config: $ex")
    }
  }

  it should "parse correct missing reference-disk-localization-manifest-files config" in {
    val backendConfig = ConfigFactory.parseString(configString())

    val validatedGcsPathsToReferenceDiskManifestFilesErrorOr = PipelinesApiConfigurationAttributes.validateGcsPathsToReferenceDiskManifestFiles(backendConfig)
    validatedGcsPathsToReferenceDiskManifestFilesErrorOr match {
      case Valid(validatedGcsPathsToReferenceDiskManifestFilesOpt) =>
        validatedGcsPathsToReferenceDiskManifestFilesOpt shouldBe None
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to reference disk manifest files from config: $ex")
    }
  }

  it should "parse correct empty reference-disk-localization-manifest-files config" in {
    val referenceDiskManifestConfigStr = "reference-disk-localization-manifest-files = []"
    val backendConfig = ConfigFactory.parseString(configString(referenceDiskManifestConfigStr))

    val validatedGcsPathsToReferenceDiskManifestFilesErrorOr = PipelinesApiConfigurationAttributes.validateGcsPathsToReferenceDiskManifestFiles(backendConfig)
    validatedGcsPathsToReferenceDiskManifestFilesErrorOr match {
      case Valid(validatedGcsPathsToReferenceDiskManifestFilesOpt) =>
        validatedGcsPathsToReferenceDiskManifestFilesOpt match {
          case Some(validatedGcsPathsToReferenceDiskManifestFiles) =>
            validatedGcsPathsToReferenceDiskManifestFiles.isEmpty shouldBe true
          case None =>
            fail("GCS paths to reference disk manifest files, parsed from config, should not be None")
        }
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to reference disk manifest files from config: $ex")
    }
  }

  it should "parse correct existing docker-image-cache-manifest-file config" in {
    val dockerImageCacheManifest1Path = "gs://bucket/manifest1.json"
    val dockerImageCacheManifestConfigStr = s"""docker-image-cache-manifest-file = "$dockerImageCacheManifest1Path""""
    val backendConfig = ConfigFactory.parseString(configString(dockerImageCacheManifestConfigStr))

    val validatedGcsPathToDockerImageCacheManifestFileErrorOr = PipelinesApiConfigurationAttributes.validateGcsPathToDockerImageCacheManifestFile(backendConfig)
    validatedGcsPathToDockerImageCacheManifestFileErrorOr match {
      case Valid(validatedGcsPathToDockerImageCacheManifestFileOpt) =>
        validatedGcsPathToDockerImageCacheManifestFileOpt match {
          case Some(validatedGcsPathToDockerCacheManifestFile) =>
            validatedGcsPathToDockerCacheManifestFile shouldBe GcsPathBuilder.validateGcsPath(dockerImageCacheManifest1Path)
          case None =>
            fail("GCS paths to docker image cache manifest files, parsed from config, should not be empty")
        }
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to docker image cache manifest files from config: $ex")
    }
  }

  it should "parse correct missing docker-image-cache-manifest-file config" in {
    val backendConfig = ConfigFactory.parseString(configString())

    val validatedGcsPathsToDockerImageCacheManifestFilesErrorOr = PipelinesApiConfigurationAttributes.validateGcsPathsToReferenceDiskManifestFiles(backendConfig)
    validatedGcsPathsToDockerImageCacheManifestFilesErrorOr match {
      case Valid(validatedGcsPathsToDockerImageCacheManifestFilesOpt) =>
        validatedGcsPathsToDockerImageCacheManifestFilesOpt shouldBe None
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to docker image cache manifest files from config: $ex")
    }
  }
}
