package cromwell.backend.google.batch.models

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigFactory}
import common.assertion.CromwellTimeoutSpec
import common.exception.MessageAggregation
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes._
import cromwell.backend.google.batch.models.GcpBatchTestConfig.BatchGlobalConfig
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.MockAuthMode
import cromwell.filesystems.gcs.GcsPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

//import java.net.URL
import scala.concurrent.duration._

class GcpBatchConfigurationAttributesSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks {

  behavior of "GcpBatchAttributes"

  val googleConfig: GoogleConfiguration = GoogleConfiguration(BatchGlobalConfig)
  val runtimeConfig: Config = ConfigFactory.load()

  it should "parse correct Batch config" in {

    val backendConfig = ConfigFactory.parseString(configString())
    println(backendConfig)

    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
    println(gcpBatchAttributes)
    gcpBatchAttributes.project should be("myProject")
    gcpBatchAttributes.executionBucket should be("gs://myBucket")
    gcpBatchAttributes.maxPollingInterval should be(600)
    gcpBatchAttributes.computeServiceAccount should be("default")
    gcpBatchAttributes.restrictMetadataAccess should be(false)
    gcpBatchAttributes.referenceFileToDiskImageMappingOpt.isEmpty should be(true)
  }

  it should "parse correct preemptible config" in {

    val backendConfig = ConfigFactory.parseString(configString(customContent = "preemptible = 3"))

    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
    gcpBatchAttributes.project should be("myProject")
    gcpBatchAttributes.executionBucket should be("gs://myBucket")
    gcpBatchAttributes.maxPollingInterval should be(600)
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
    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")

    gcpBatchAttributes.batchRequestTimeoutConfiguration.readTimeoutMillis.get.value should be(100.hours.toMillis.toInt)
    gcpBatchAttributes.batchRequestTimeoutConfiguration.connectTimeoutMillis.get.value should be(
      10.seconds.toMillis.toInt
    )
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
    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")

    gcpBatchAttributes.batchRequestTimeoutConfiguration should be(BatchRequestTimeoutConfiguration(None, None))
  }

  it should "parse batch-timeout" in {

    val backendConfig = ConfigFactory.parseString(configString(customContent = "batch-timeout = 3 days"))
    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")

    gcpBatchAttributes.batchTimeout should be(3.days)
  }

  it should "parse an undefined batch-timeout" in {

    val backendConfig = ConfigFactory.parseString(configString())
    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")

    gcpBatchAttributes.batchTimeout should be(7.days)
  }

  it should "parse compute service account" in {

    val backendConfig = ConfigFactory.parseString(configString(genomics = """compute-service-account = "testing" """))

    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
    gcpBatchAttributes.computeServiceAccount should be("testing")
  }

  it should "parse restrict-metadata-access" in {

    val backendConfig = ConfigFactory.parseString(configString(genomics = "restrict-metadata-access = true"))

    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
    gcpBatchAttributes.restrictMetadataAccess should be(true)

  }

  it should "parse localization-attempts" in {

    val backendConfig = ConfigFactory.parseString(configString(genomics = "localization-attempts = 31380"))

    val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
    gcpBatchAttributes.gcsTransferConfiguration.transferAttempts.value should be(31380)
  }

  private val mockAuth = MockAuthMode("mock")

  private val validVpcConfigTests = Table(
    ("description", "customConfig", "vpcConfig"),
    ("empty body", "virtual-private-cloud {}", VirtualPrivateCloudConfiguration(None, None)),
    (
      "labels config",
      """virtual-private-cloud {
        |  network-label-key = my-network
        |  subnetwork-label-key = my-subnetwork
        |  auth = mock
        |}
        |""".stripMargin,
      VirtualPrivateCloudConfiguration(
        Option(VirtualPrivateCloudLabels("my-network", Option("my-subnetwork"), mockAuth)),
        None
      )
    ),
    (
      "labels config without subnetwork key",
      """virtual-private-cloud {
        |  network-label-key = my-network
        |  auth = mock
        |}
        |""".stripMargin,
      VirtualPrivateCloudConfiguration(
        Option(VirtualPrivateCloudLabels("my-network", None, mockAuth)),
        None
      )
    ),
    (
      "literal config",
      """virtual-private-cloud {
        |  network-name = my-network
        |  subnetwork-name = my-subnetwork
        |}
        |""".stripMargin,
      VirtualPrivateCloudConfiguration(
        None,
        Option(VirtualPrivateCloudLiterals("my-network", Option("my-subnetwork")))
      )
    ),
    (
      "literal config without subnetwork name",
      """virtual-private-cloud {
        |  network-name = my-network
        |}
        |""".stripMargin,
      VirtualPrivateCloudConfiguration(
        None,
        Option(VirtualPrivateCloudLiterals("my-network", None))
      )
    )
  )

  private val invalidVPCConfigTests = Table(
    ("description", "customConfig", "messages"),
    (
      "without auth",
      """virtual-private-cloud {
        |  network-label-key = my-network
        |}
        |""".stripMargin,
      List("Virtual Private Cloud configuration is invalid. Missing keys: `auth`.")
    ),
    (
      "without network label-key",
      """virtual-private-cloud {
        |  auth = mock
        |}
        |""".stripMargin,
      List("Virtual Private Cloud configuration is invalid. Missing keys: `network-label-key`.")
    ),
    (
      "with just a subnetwork label key",
      """virtual-private-cloud {
        |  subnetwork-label-key = my-subnetwork
        |}
        |""".stripMargin,
      List("Virtual Private Cloud configuration is invalid. Missing keys: `network-label-key,auth`.")
    ),
    (
      "with subnetwork label network key and auth",
      """ virtual-private-cloud {
        |   subnetwork-label-key = my-subnetwork
        |   auth = mock
        | }
        |""".stripMargin,
      List("Virtual Private Cloud configuration is invalid. Missing keys: `network-label-key`.")
    )
  )

  forAll(validVpcConfigTests) { (description, customConfig, vpcConfig) =>
    it should s"parse virtual-private-cloud $description" in {

      val backendConfig = ConfigFactory.parseString(configString(customConfig))
      val gcpBatchAttributes = GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
      gcpBatchAttributes.virtualPrivateCloudConfiguration should be(vpcConfig)
    }
  }

  forAll(invalidVPCConfigTests) { (description, customConfig, errorMessages) =>
    it should s"not parse invalid virtual-private-cloud config $description" in {

      val backendConfig = ConfigFactory.parseString(configString(customConfig))
      val exception = intercept[IllegalArgumentException with MessageAggregation] {
        GcpBatchConfigurationAttributes(googleConfig, backendConfig, "batch")
      }
      exception.errorMessages.toList should be(errorMessages)
    }
  }

  it should "not parse invalid config" in {

    val nakedConfig =
      ConfigFactory.parseString("""
                                  |{
                                  |   genomics {
                                  |
                                  |   }
                                  |}
        """.stripMargin)

    val exception = intercept[IllegalArgumentException with MessageAggregation] {
      GcpBatchConfigurationAttributes(googleConfig, nakedConfig, "batch")
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("String: 2: No configuration setting found for key 'project'")
    errorsList should contain("String: 2: No configuration setting found for key 'root'")
    errorsList should contain("String: 3: No configuration setting found for key 'genomics.auth'")
    errorsList should contain("String: 2: No configuration setting found for key 'filesystems'")
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
       |     auth = "mock"
       |    $genomics
       |     endpoint-url = "http://myEndpoint"
       |     location = "us-central1"
       |   }
       |
       |   filesystems = {
       |     gcs {
       |       // A reference to a potentially different auth for manipulating files via engine functions.
       |       auth = "mock"
       |     }
       |   }
       |}
       | """.stripMargin

  it should "parse gsutil memory specifications" in {

    val valids = List("0", "150M", "14   PIBIT", "6kib")

    valids foreach {
      case GcpBatchConfigurationAttributes.GsutilHumanBytes(_, _) =>
      case bad => fail(s"'$bad' was expected to be a valid gsutil memory specification")
    }
  }

  it should "reject invalid memory specifications" in {

    val invalids = List("-1", "150MB", "14PB")

    invalids foreach {
      case invalid @ GcpBatchConfigurationAttributes.GsutilHumanBytes(_, _) =>
        fail(s"Memory specification $invalid not expected to be accepted")
      case _ =>
    }
  }

  it should "parse a missing \"reference-disk-localization-manifests\"" in {

    val backendConfig = ConfigFactory.parseString(configString())

    val validation = GcpBatchConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "batch")

    validation shouldBe None.validNel
  }

  it should "parse a present but empty \"reference-disk-localization-manifests\"" in {

    val manifestConfig = "reference-disk-localization-manifests = []"

    val backendConfig = ConfigFactory.parseString(configString(customContent = manifestConfig))

    val validation = GcpBatchConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "batch")

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
    val validation = GcpBatchConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "batch")
    val manifests: List[ManifestFile] = validation.toEither.toOption.get.get

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
         |}]"""
    )

    badValues foreach { badValue =>
      val customContent = s""""reference-disk-localization-manifests" = $badValue"""
      val backendConfig = ConfigFactory.parseString(configString(customContent))
      val validation = GcpBatchConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "batch")
      validation.isInvalid shouldBe true
    }
  }

  it should "parse correct existing docker-image-cache-manifest-file config" in {

    val dockerImageCacheManifest1Path = "gs://bucket/manifest1.json"
    val dockerImageCacheManifestConfigStr = s"""docker-image-cache-manifest-file = "$dockerImageCacheManifest1Path""""
    val backendConfig = ConfigFactory.parseString(configString(dockerImageCacheManifestConfigStr))

    val validatedGcsPathToDockerImageCacheManifestFileErrorOr =
      GcpBatchConfigurationAttributes.validateGcsPathToDockerImageCacheManifestFile(backendConfig)
    validatedGcsPathToDockerImageCacheManifestFileErrorOr match {
      case Valid(validatedGcsPathToDockerImageCacheManifestFileOpt) =>
        validatedGcsPathToDockerImageCacheManifestFileOpt match {
          case Some(validatedGcsPathToDockerCacheManifestFile) =>
            validatedGcsPathToDockerCacheManifestFile shouldBe GcsPathBuilder.validateGcsPath(
              dockerImageCacheManifest1Path
            )
          case None =>
            fail("GCS paths to docker image cache manifest files, parsed from config, should not be empty")
        }
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to docker image cache manifest files from config: $ex")
    }
  }

  it should "parse correct missing docker-image-cache-manifest-file config" in {

    val backendConfig = ConfigFactory.parseString(configString())

    val validatedGcsPathsToDockerImageCacheManifestFilesErrorOr =
      GcpBatchConfigurationAttributes.validateReferenceDiskManifestConfigs(backendConfig, "unit-test-backend")
    validatedGcsPathsToDockerImageCacheManifestFilesErrorOr match {
      case Valid(validatedGcsPathsToDockerImageCacheManifestFilesOpt) =>
        validatedGcsPathsToDockerImageCacheManifestFilesOpt shouldBe None
      case Invalid(ex) =>
        fail(s"Error while parsing GCS paths to docker image cache manifest files from config: $ex")
    }
  }
}
