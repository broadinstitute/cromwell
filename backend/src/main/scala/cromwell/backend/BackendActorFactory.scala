package cromwell.backend

import akka.actor.{ActorRef, ActorSystem}
import cromwell.backend.model.WorkflowDescriptor

/**
  * Creates Backend instances based on backend name.
  */
trait BackendActorFactory {
  /**
    * Returns a backend actor instance based on the name, Cromwell actor system and specific task.
    *
    * @param initClass Initialization class for the specific backend.
    * @param actorSystem Cromwell Engine actor system.
    * @param workflowDescriptor Needed data to be able to execute a workflow.
    * @return An actor of type BackendActor.
    */
  def getBackend(initClass: String, actorSystem: ActorSystem, workflowDescriptor: WorkflowDescriptor): ActorRef
}