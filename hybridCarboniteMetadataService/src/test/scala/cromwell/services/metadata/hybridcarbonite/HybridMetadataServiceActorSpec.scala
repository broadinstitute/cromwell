package cromwell.services.metadata.hybridcarbonite

import akka.testkit.TestProbe
import cromwell.services.metadata.impl.MetadataServiceActorSpec
import cromwell.services.metadata.impl.MetadataServiceActorSpec.globalConfigToMetadataServiceConfig

class HybridMetadataServiceActorSpec extends MetadataServiceActorSpec {

  override def actorName: String = "HybridMetadataServiceActor"
  override val actor = system.actorOf(HybridMetadataServiceActor.props(config, globalConfigToMetadataServiceConfig(config), TestProbe().ref), "HybridMetadataServiceActor-for-MetadataServiceActorSpec")

}
