package cromwell.services

import com.typesafe.config.ConfigFactory
import cromwell.server.WorkflowManagerSystem

private object ServiceRegistryInstance extends WorkflowManagerSystem {
  val ServiceRegistryActorInstance = actorSystem.actorOf(ServiceRegistryActor.props(ConfigFactory.load()), "ServiceRegistryActor")
}

trait ServiceRegistryClient {
  def serviceRegistryActor = ServiceRegistryInstance.ServiceRegistryActorInstance
}
