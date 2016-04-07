package cromwell.backend

import akka.actor.{ActorRef, ActorSystem}
import wdl4s.Call

trait WorkflowBackendActorFactory {

  /**
    * Creates a workflow backend actor ref based on the init class, Cromwell actor system,
    * backend workflow descriptor, backend configuration and Cromwell service registry.
    *
    * @param initClass                      Initialization class for the specific backend.
    * @param actorSystem                    Cromwell Engine actor system.
    * @param backendWorkflowDescriptor      Workflow level data needed to be able to execute a workflow.
    * @param backendConfigurationDescriptor Backend configuration data.
    * @param serviceRegistryActor           Cromwell service registry.
    * @param calls                          A sequence of WDL Calls used to validate runtime requirements by the backend.
    * @return An actor ref pointing to a WorkflowBackendActor implementation.
    */
  def getBackend(initClass: String,
                 actorSystem: ActorSystem,
                 backendWorkflowDescriptor: BackendWorkflowDescriptor,
                 backendConfigurationDescriptor: BackendConfigurationDescriptor,
                 serviceRegistryActor: ActorRef,
                 calls: Seq[Call]): ActorRef

}
