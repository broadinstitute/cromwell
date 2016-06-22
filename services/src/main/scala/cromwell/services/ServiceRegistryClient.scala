package cromwell.services

import com.typesafe.config.ConfigFactory

private object ServiceRegistryInstance {
  // TODO: PBE: Removed engine's WMS trait creating a singleton actor system... can't we just pass actor refs?
  println("CLASSPATH is " + getClass.getClassLoader)
  private val actorSystem = akka.actor.ActorSystem("cromwell-service-registry-system")
  val ServiceRegistryActorInstance = actorSystem.actorOf(ServiceRegistryActor.props(ConfigFactory.load()), "ServiceRegistryActor")
}

trait ServiceRegistryClient {
  // This is a `val` rather than the expected `def` to force the ServiceRegistryActor to initialize.  This catches
  // config errors up front and updates the DB schema promptly after startup.  If this is a `def` or `lazy val` these
  // things only happen the first time an endpoint is accessed.
  val serviceRegistryActor = ServiceRegistryInstance.ServiceRegistryActorInstance
}
