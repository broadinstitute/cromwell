package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.standard._
import wdl4s.TaskCall

case class TesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override def initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[TesInitializationActor]

  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[TesAsyncBackendJobExecutionActor]

  override def jobIdKey: String = TesAsyncBackendJobExecutionActor.JobIdKey

  val tesConfiguration = new TesConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                                 serviceRegistryActor: ActorRef): StandardInitializationActorParams = {
    TesInitializationActorParams(workflowDescriptor, calls, tesConfiguration, serviceRegistryActor)
  }
}
