package cromwell.backend.impl.bcs

import akka.actor.ActorRef
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.backend.standard._
import wom.graph.CommandCallNode

final case class BcsBackendLifecycleActorFactory(val name: String, val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {
  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[BcsInitializationActor]
  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] = classOf[BcsAsyncBackendJobExecutionActor]

  override def jobIdKey: String = BcsAsyncBackendJobExecutionActor.JobIdKey

  val bcsConfiguration = new BcsConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode], serviceRegistryActor: ActorRef, restarting: Boolean): StandardInitializationActorParams = {
    BcsInitializationActorParams(workflowDescriptor, calls, bcsConfiguration, serviceRegistryActor)
  }
}
