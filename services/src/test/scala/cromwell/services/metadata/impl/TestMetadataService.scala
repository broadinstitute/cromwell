package cromwell.services.metadata.impl

import akka.actor.ActorRef
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.services.ServiceRegistryActor

object TestMetadataService {
  private lazy val ServiceRegistryActorSystem = akka.actor.ActorSystem("cromwell-test-metadata-summary-system")

  lazy val ServiceRegistryActorInstance: ActorRef = {
    ServiceRegistryActorSystem.actorOf(ServiceRegistryActor.props(ConfigFactory.load()), "TestMetadataSummaryActor")
  }
}

class TestMetadataService(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends MetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor) {
  override private [impl] def buildSummaryActor = Option(TestMetadataService.ServiceRegistryActorInstance)
}
