package cromwell.backend.impl.tes

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.core.Dispatcher._
import wdl4s.TaskCall

case class TesBackendLifecycleActorFactory(name: String,
                                           configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Set[TaskCall],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(TesInitializationActor.props(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor).withDispatcher(BackendDispatcher))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    TesJobExecutionActor.props(jobDescriptor, configurationDescriptor).withDispatcher(BackendDispatcher)
  }
}
