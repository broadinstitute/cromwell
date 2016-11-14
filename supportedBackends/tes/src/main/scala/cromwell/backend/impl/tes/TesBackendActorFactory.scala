package cromwell.backend.impl.tes

import akka.actor.{ActorRef, Props}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData,
                         BackendJobDescriptor, BackendLifecycleActorFactory,
                         BackendWorkflowDescriptor}
import wdl4s.Call

case class TesBackendActorFactory(name: String,
                                  configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(TesInitializationActor.props(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    TesJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }
}
