package cromwell.services

import com.typesafe.config.ConfigFactory

private object ServiceRegistryInstance {
  // TODO: PBE: Removed engine's WMS trait creating a singleton actor system... can't we just pass actor refs?
  private val actorSystem = akka.actor.ActorSystem("cromwell-service-registry-system")
  val ServiceRegistryActorInstance = actorSystem.actorOf(ServiceRegistryActor.props(ConfigFactory.load()), "ServiceRegistryActor")
}

trait ServiceRegistryClient {
  def serviceRegistryActor = ServiceRegistryInstance.ServiceRegistryActorInstance
}
