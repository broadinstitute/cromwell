package cromwell.backend.sfs.callcaching

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.routing.RouterConfig
import cromwell.backend.callcaching.BackendHashingMethods
import cromwell.backend.validation.RuntimeAttributesKeys

case class ConfigSfsBackendHashingMethods(actorSystem: ActorSystem) extends BackendHashingMethods {

  // Docker not required since it's already handed by the engine. Hooray!
  override def hashableRuntimeAttributes: List[String] = List(
    RuntimeAttributesKeys.ContinueOnReturnCodeKey,
    RuntimeAttributesKeys.CpuKey,
    RuntimeAttributesKeys.FailOnStderrKey,
    RuntimeAttributesKeys.MemoryKey
  )

  override val fileContentsHasherActor: ActorRef = actorSystem.actorOf(Props(new ConfigSfsFileContentsHasherActor))
}
