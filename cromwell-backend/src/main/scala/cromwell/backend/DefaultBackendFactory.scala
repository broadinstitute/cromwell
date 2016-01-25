package cromwell.backend

import akka.actor.{ActorRef, ActorSystem, Props}
import cromwell.backend.model.TaskDescriptor

object DefaultBackendFactory extends BackendFactory {
  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    * @param initClass Initialization class for the specific backend.
    * @param actorSystem Cromwell Engine actor system.
    * @param task Specific task to be executed in the backend.
    * @return A backend instance.
    */
  override def getBackend(initClass: String, actorSystem: ActorSystem, task: TaskDescriptor): ActorRef = {
    val backendActor = Props.create(Class.forName(initClass).asInstanceOf[Class[BackendActor]], task)
    actorSystem.actorOf(backendActor)
  }
}
