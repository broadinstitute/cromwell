package cromwell.backend

import akka.actor.{ActorRef, ActorSystem}

/**
  * Creates Backend instances based on backend name.
  */
trait BackendActorFactory {
  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    *
    * @param initClass Initialization class for the specific backend.
    * @param actorSystem Cromwell Engine actor system.
    * @return An actor of type BackendActor.
    */
  def getBackend(initClass: String, actorSystem: ActorSystem): ActorRef
}