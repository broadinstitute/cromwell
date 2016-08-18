package cromwell.backend.impl.jes.callcaching

import akka.actor.{ActorRef, ActorSystem, Props}
import cromwell.backend.callcaching.BackendHashingMethods
import cromwell.backend.impl.jes.JesRuntimeAttributes
import cromwell.backend.validation.RuntimeAttributesKeys
import cromwell.core.callcaching.HashKey

case class JesBackendHashingMethods(actorSystem: ActorSystem) extends BackendHashingMethods {

  // Docker not required since it's already handed by the engine. Hooray!
  override def hashableRuntimeAttributes: List[String] = List(
    RuntimeAttributesKeys.ContinueOnReturnCodeKey,
    RuntimeAttributesKeys.CpuKey,
    JesRuntimeAttributes.DisksKey,
    JesRuntimeAttributes.ZonesKey,
    RuntimeAttributesKeys.FailOnStderrKey,
    RuntimeAttributesKeys.MemoryKey,
    JesRuntimeAttributes.PreemptibleKey,
    JesRuntimeAttributes.BootDiskSizeKey
  )

  override val fileContentsHasherActor: ActorRef = actorSystem.actorOf(Props(new JesFileContentsHasherActor))
}
