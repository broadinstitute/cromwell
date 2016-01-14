package cromwell.engine.backend

import akka.actor.ActorSystem
import cromwell.backend.{Backend, BackendActor, BackendFactory}
import cromwell.backend.model.TaskDescriptor

class DefaultBackendFactory extends BackendFactory {
  /**
    * Returns a backend instance based on the name and specific task.
    * @param name Backend name.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  override def getBackend(name: String, task: TaskDescriptor): cromwell.backend.Backend = ???

  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    * @param name Backend name.
    * @param actorSystem Cromwell Engine actor system.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  override def getBackend(name: String, actorSystem: ActorSystem, task: TaskDescriptor): BackendActor = ???
}
