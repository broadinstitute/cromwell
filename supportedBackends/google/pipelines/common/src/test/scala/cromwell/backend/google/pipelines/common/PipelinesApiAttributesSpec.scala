package cromwell.backend.google.pipelines.common

import java.net.URL

import com.typesafe.config.ConfigFactory
import common.exception.MessageAggregation
import cromwell.cloudsupport.gcp.GoogleConfiguration
import org.scalatest.{FlatSpec, Matchers}

class PipelinesApiAttributesSpec extends FlatSpec with Matchers {

  import PipelinesApiTestConfig._

  behavior of "PipelinesApiAttributes"

  val googleConfig = GoogleConfiguration(PapiGlobalConfig)
  val runtimeConfig = ConfigFactory.load()

  it should "parse correct PAPI config" in {

    val backendConfig = ConfigFactory.parseString(configString())

    val pipelinesApiAttributes = PipelinesApiAttributes(googleConfig, backendConfig)
    pipelinesApiAttributes.endpointUrl should be(new URL("http://myEndpoint"))
    pipelinesApiAttributes.project should be("myProject")
    pipelinesApiAttributes.executionBucket should be("gs://myBucket")
    pipelinesApiAttributes.maxPollingInterval should be(600)
    pipelinesApiAttributes.computeServiceAccount should be("default")
    pipelinesApiAttributes.restrictMetadataAccess should be(false)
  }

  it should "parse correct preemptible config" in {
    val backendConfig = ConfigFactory.parseString(configString(preemptible = "preemptible = 3"))

    val pipelinesApiAttributes = PipelinesApiAttributes(googleConfig, backendConfig)
    pipelinesApiAttributes.endpointUrl should be(new URL("http://myEndpoint"))
    pipelinesApiAttributes.project should be("myProject")
    pipelinesApiAttributes.executionBucket should be("gs://myBucket")
    pipelinesApiAttributes.maxPollingInterval should be(600)
  }

  it should "parse compute service account" in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = """compute-service-account = "testing" """))

    val pipelinesApiAttributes = PipelinesApiAttributes(googleConfig, backendConfig)
    pipelinesApiAttributes.computeServiceAccount should be("testing")
  }

  it should "parse restrict-metadata-access" in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = "restrict-metadata-access = true"))

    val pipelinesApiAttributes = PipelinesApiAttributes(googleConfig, backendConfig)
    pipelinesApiAttributes.restrictMetadataAccess should be(true)

  }

  it should "parse localization-attempts" in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = "localization-attempts = 31380"))

    val pipelinesApiAttributes = PipelinesApiAttributes(googleConfig, backendConfig)
    pipelinesApiAttributes.localizationConfiguration.localizationAttempts.value should be(31380)

  }

  it should "parse virtual-private-cloud" in {
    val backendConfig = ConfigFactory.parseString(configString(networkLabelKey = "network-label-key = my-network", vpcAuth = "auth = application-default"))

    val pipelinesApiAttributes = PipelinesApiAttributes(googleConfig, backendConfig)
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.name should be("my-network")
    pipelinesApiAttributes.virtualPrivateCloudConfiguration.get.auth.name should be("application-default")
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
      PipelinesApiAttributes(googleConfig, nakedConfig)
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("No configuration setting found for key 'project'")
    errorsList should contain("No configuration setting found for key 'root'")
    errorsList should contain("No configuration setting found for key 'genomics.auth'")
    errorsList should contain("No configuration setting found for key 'filesystems'")
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
      PipelinesApiAttributes(googleConfig, networkKeyOnlyConfig)
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Auth scheme not provided for Virtual Private Cloud configuration.")
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
      PipelinesApiAttributes(googleConfig, authOnlyConfig)
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("Network label key not provided for Virtual Private Cloud configuration.")
  }

  def configString(preemptible: String = "", genomics: String = "", networkLabelKey: String = "", vpcAuth: String = ""): String =
    s"""
      |{
      |   project = "myProject"
      |   root = "gs://myBucket"
      |   maximum-polling-interval = 600
      |
      |   virtual-private-cloud {
      |     $networkLabelKey
      |     $vpcAuth
      |   }
      |
      |   $preemptible
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
}
