package cromwell.services.metadata.hybridcarbonite

import akka.testkit.TestProbe
import com.typesafe.config.ConfigFactory
import cromwell.services.metadata.impl.MetadataServiceActorSpec
import cromwell.services.metadata.impl.MetadataServiceActorSpec.globalConfigToMetadataServiceConfig

class HybridMetadataServiceActorSpec extends MetadataServiceActorSpec {

  override def actorName: String = "HybridMetadataServiceActor"

  val hybridConfigString =
    s"""${HybridMetadataServiceActor.CarboniteConfigPath} {
        |    carbonite-interval = Inf
        |
        |    bucket = "this test shouldn't need a bucket"
        |
        |    filesystems {
        |      # This wouldn't work in the real world... this is just for the tests
        |      local {
        |      }
        |    }
        |}""".stripMargin

  val hybridConfig = ConfigFactory.parseString(hybridConfigString)

  override lazy val actor = system.actorOf(HybridMetadataServiceActor.props(hybridConfig, globalConfigToMetadataServiceConfig(config), TestProbe().ref), "HybridMetadataServiceActor-for-MetadataServiceActorSpec")

}
