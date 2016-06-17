package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.filesystems.gcs.GoogleConfiguration
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.ThrowableWithErrors

class JesAttributesSpec extends FlatSpec with Matchers {

  it should "parse correct JES config" in {
//    val configString = """
//          {
//             project = "myProject"
//             root = "gs://myBucket"
//             maximum-polling-interval = 600
//             [PREEMPTIBLE]
//             genomics {
//               // A reference to an auth defined in the `google` stanza at the top.  This auth is used to create
//               // Pipelines and manipulate auth JSONs.
//               auth = "service-account"
//               endpoint-url = "http://myEndpoint"
//             }
//
//             filesystems = {
//               gcs {
//                 // A reference to a potentially different auth for manipulating files via engine functions.
//                 auth = "service-account"
//               }
//             }
//          }""".stripMargin
//
//    val fullConfig = ConfigFactory.parseString(configString.replace("[PREEMPTIBLE]", "preemptible = 3"))
//    val googleConfig = GoogleConfiguration(ConfigFactory.load)
//
//    val fullAttributes = JesAttributes.apply(googleConfig, BackendConfigurationDescriptor(fullConfig, ConfigFactory.load()))
//    fullAttributes.endpointUrl.toString shouldBe new URL("http://myEndpoint").toString
//    fullAttributes.project shouldBe "myProject"
//    fullAttributes.executionBucket shouldBe "gs://myBucket"
//    fullAttributes.maxPollingInterval shouldBe 600
//
//    val noPreemptibleConfig = ConfigFactory.parseString(configString.replace("[PREEMPTIBLE]", ""))
//
//    val noPreemptibleAttributes = JesAttributes.apply(googleConfig, BackendConfigurationDescriptor(noPreemptibleConfig, ConfigFactory.load))
//    noPreemptibleAttributes.endpointUrl.toString shouldBe new URL("http://myEndpoint").toString
//    noPreemptibleAttributes.project shouldBe "myProject"
//    noPreemptibleAttributes.executionBucket shouldBe "gs://myBucket"
//    noPreemptibleAttributes.maxPollingInterval shouldBe 600
//  }
//
//  it should "not parse invalid config" in {
//    val nakedConfig =
//      ConfigFactory.parseString("""
//        |{
//        |   genomics {
//        |     endpoint-url = "myEndpoint"
//        |   }
//        |}
//      """.stripMargin)
//
//    val googleConfig = GoogleConfiguration(ConfigFactory.load)
//
//    val exception = intercept[IllegalArgumentException with ThrowableWithErrors] {
//      JesAttributes.apply(googleConfig, BackendConfigurationDescriptor(nakedConfig, ConfigFactory.load))
//    }
//    val errorsList = exception.errors.list
//    errorsList should contain ("Could not find key: project")
//    errorsList should contain ("Could not find key: root")
//    errorsList should contain ("no protocol: myEndpoint")
  }

}
