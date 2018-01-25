package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.standard._
import wom.graph.CommandCallNode

case class TesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[TesInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[TesAsyncBackendJobExecutionActor]

  override def jobIdKey: String = TesAsyncBackendJobExecutionActor.JobIdKey

  val tesConfiguration = new TesConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                                 serviceRegistryActor: ActorRef, restarting: Boolean): StandardInitializationActorParams = {
    TesInitializationActorParams(workflowDescriptor, calls, tesConfiguration, serviceRegistryActor)
  }
}
