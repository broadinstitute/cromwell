package cromwell.backend

import akka.actor.{ActorRef, ActorSystem}
import cromwell.backend.model.TaskDescriptor

/**
  * Creates Backend instances based on backend name, actor system and task.
  */
trait BackendFactory {
  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    * @param initClass Initialization class for the specific backend.
    * @param actorSystem Cromwell Engine actor system.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  def getBackend(initClass: String, actorSystem: ActorSystem, task: TaskDescriptor): ActorRef
}
