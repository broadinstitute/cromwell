package cromwell.backend.config

import akka.actor.ActorSystem
import cromwell.backend.model.TaskDescriptor
import cromwell.backend.{Backend, BackendActor, BackendFactory}

class DefaultBackendFactory extends BackendFactory {
  /**
    * Returns a backend instance based on the name and specific task.
    * @param initClass Initialization class for the specific backend.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  override def getBackend(initClass: String, task: TaskDescriptor): Backend = {
    val constructor = Class.forName(initClass).getConstructor(TaskDescriptor.getClass)
    constructor.newInstance(task).asInstanceOf[Backend]
  }

  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    * @param initClass Initialization class for the specific backend.
    * @param actorSystem Cromwell Engine actor system.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  override def getBackend(initClass: String, actorSystem: ActorSystem, task: TaskDescriptor): BackendActor = {
    val constructor = Class.forName(initClass).getConstructor(ActorSystem.getClass, TaskDescriptor.getClass)
    constructor.newInstance(actorSystem, task).asInstanceOf[BackendActor]
  }
}
