package cromwell.backend.google.pipelines.common

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import common.exception.MessageAggregation
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.BatchRequestTimeoutConfiguration
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.filesystems.gcs.GcsPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.net.URL
import scala.concurrent.duration._

class PipelinesApiConfigurationAttributesSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  import PipelinesApiTestConfig._

  behavior of "PipelinesApiAttributes"

  val googleConfig: GoogleConfiguration = GoogleConfiguration(PapiGlobalConfig)
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

  it should "parse a missing \"reference-disk-localization-manifests\"" in {
    val backendConfig = ConfigFactory.parseString(configString())

    val validation = PipelinesApiConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "papi")

    validation shouldBe None.validNel
  }

  it should "parse a present but empty \"reference-disk-localization-manifests\"" in {
    val manifestConfig = "reference-disk-localization-manifests = []"

    val backendConfig = ConfigFactory.parseString(configString(customContent = manifestConfig))

    val validation = PipelinesApiConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "papi")

    validation shouldBe Option(List.empty).validNel
  }

  it should "parse a present and populated \"reference-disk-localization-manifests\"" in {
    // Highly abridged versions of hg19 and hg38 manifests just to test for correctness
    // of parsing.
    val manifestConfig =
    """
      |reference-disk-localization-manifests = [
      |{
      |  "imageIdentifier" : "hg19-public-2020-10-26",
      |  "diskSizeGb" : 10,
      |  "files" : [ {
      |    "path" : "gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
      |    "crc32c" : 159565724
      |  }, {
      |    "path" : "gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict",
      |    "crc32c" : 1679459712
      |  }]
      |},
      |{
      |  "imageIdentifier" : "hg38-public-2020-10-26",
      |  "diskSizeGb" : 20,
      |  "files" : [ {
      |    "path" : "gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
      |    "crc32c" : 930173616
      |  }, {
      |    "path" : "gcp-public-data--broad-references/hg38/v0/exome_evaluation_regions.v1.interval_list",
      |    "crc32c" : 289077232
      |  }]
      |}
      |]
      |""".stripMargin
    val backendConfig = ConfigFactory.parseString(configString(manifestConfig))
    val validation = PipelinesApiConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "papi")
    val manifests: List[ManifestFile] = validation.toEither.right.get.get

    manifests shouldBe List(
      ManifestFile(
        imageIdentifier = "hg19-public-2020-10-26",
        diskSizeGb = 10,
        files = List(
          ReferenceFile(
            path = "gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
            crc32c = 159565724
          ),
          ReferenceFile(
            path = "gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict",
            crc32c = 1679459712
          )
        )
      ),
      ManifestFile(
        imageIdentifier = "hg38-public-2020-10-26",
        diskSizeGb = 20,
        files = List(
          ReferenceFile(
            path = "gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            crc32c = 930173616
          ),
          ReferenceFile(
            path = "gcp-public-data--broad-references/hg38/v0/exome_evaluation_regions.v1.interval_list",
            crc32c = 289077232
          )
        )
      )
    )
  }

  it should "parse a present and invalid \"reference-disk-localization-manifests\"" in {
    val badValues = List(
      "\"foo\"",
      "{ foo: bar }",
      s"""
         |[{
         |   # missing imageIdentifier
         |   "diskSizeGb" : 10,
         |   "files" : [ {
         |     "path" : "gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
         |     "crc32c" : 159565724
         |     }]
         |}]""",
      s"""
         |[{
         |   "imageIdentifier" : "hg19-public-2020-10-26",
         |   # missing diskSizeGb
         |   "files" : [ {
         |     "path" : "gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
         |     "crc32c" : 159565724
         |     }]
         |}]""",
      s"""
         |[{
         |   "imageIdentifier" : "hg19-public-2020-10-26",
         |   "diskSizeGb" : 10,
         |   # missing files
         |}]""",
    )

    badValues foreach { badValue =>
      val customContent = s""""reference-disk-localization-manifests" = $badValue"""
      val backendConfig = ConfigFactory.parseString(configString(customContent))
      val validation = PipelinesApiConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "papi")
      validation.isInvalid shouldBe true
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

    val validatedGcsPathsToDockerImageCacheManifestFilesErrorOr = PipelinesApiConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "unit-test-backend")
    validatedGcsPathsToDockerImageCacheManifestFilesErrorOr match {
      case Valid(validatedGcsPathsToDockerImageCacheManifestFilesOpt) =>
        validatedGcsPathsToDockerImageCacheManifestFilesOpt shouldBe None
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to docker image cache manifest files from config: $ex")
    }
  }
}
