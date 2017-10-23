package cromwell.backend.impl.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.Tags._
import lenthall.exception.MessageAggregation
import org.scalatest.{FlatSpec, Matchers}

class JesAttributesSpec extends FlatSpec with Matchers {

  import JesTestConfig._

  behavior of "JesAttributes"

  val googleConfig = GoogleConfiguration(JesGlobalConfig)
  val runtimeConfig = ConfigFactory.load()

  it should "parse correct JES config" taggedAs IntegrationTest in {

    val backendConfig = ConfigFactory.parseString(configString())

    val jesAttributes = JesAttributes(googleConfig, backendConfig)
    jesAttributes.endpointUrl should be(new URL("http://myEndpoint"))
    jesAttributes.project should be("myProject")
    jesAttributes.executionBucket should be("gs://myBucket")
    jesAttributes.maxPollingInterval should be(600)
    jesAttributes.computeServiceAccount should be("default")
    jesAttributes.restrictMetadataAccess should be(false)
  }

  it should "parse correct preemptible config" taggedAs IntegrationTest in {
    val backendConfig = ConfigFactory.parseString(configString(preemptible = "preemptible = 3"))

    val jesAttributes = JesAttributes(googleConfig, backendConfig)
    jesAttributes.endpointUrl should be(new URL("http://myEndpoint"))
    jesAttributes.project should be("myProject")
    jesAttributes.executionBucket should be("gs://myBucket")
    jesAttributes.maxPollingInterval should be(600)
  }

  it should "parse compute service account" taggedAs IntegrationTest in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = """compute-service-account = "testing" """))

    val jesAttributes = JesAttributes(googleConfig, backendConfig)
    jesAttributes.computeServiceAccount should be("testing")
  }

  it should "parse restrict-metadata-access" taggedAs IntegrationTest in {
    val backendConfig = ConfigFactory.parseString(configString(genomics = "restrict-metadata-access = true"))

    val jesAttributes = JesAttributes(googleConfig, backendConfig)
    jesAttributes.restrictMetadataAccess should be(true)

  }

  it should "not parse invalid config" taggedAs IntegrationTest in {
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
      JesAttributes(googleConfig, nakedConfig)
    }
    val errorsList = exception.errorMessages.toList
    errorsList should contain("No configuration setting found for key 'project'")
    errorsList should contain("No configuration setting found for key 'root'")
    errorsList should contain("No configuration setting found for key 'genomics.auth'")
    errorsList should contain("No configuration setting found for key 'filesystems'")
    errorsList should contain("URI is not absolute")
  }

  def configString(preemptible: String = "", genomics: String = ""): String =
    s"""
      |{
      |   project = "myProject"
      |   root = "gs://myBucket"
      |   maximum-polling-interval = 600
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
