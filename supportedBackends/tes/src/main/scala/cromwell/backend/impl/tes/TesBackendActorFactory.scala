package cromwell.backend.impl.tes

import akka.actor.{ActorRef, Props}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor, BackendLifecycleActorFactory}

case class TesBackendActorFactory(name: String,
                                  configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    TesJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }
}
