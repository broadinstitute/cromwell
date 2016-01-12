package cromwell.backend

import akka.actor.ActorSystem
import cromwell.backend.model.TaskDescriptor

/**
  * Creates Backend instances based on backend name, actor system and task.
  */
trait BackendFactory {
  /**
    * Returns a backend instance based on the name and specific task.
    * @param name Backend name.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  def getBackend(name: String, task: TaskDescriptor) : Backend

  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    * @param name Backend name.
    * @param actorSystem Cromwell Engine actor system.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  def getBackend(name: String, actorSystem: ActorSystem, task: TaskDescriptor) : BackendActor
}
