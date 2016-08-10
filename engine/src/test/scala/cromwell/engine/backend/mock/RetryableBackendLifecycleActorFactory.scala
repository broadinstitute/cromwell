package cromwell.engine.backend.mock

import akka.actor.{ActorRef, ActorSystem, Props}
import cromwell.backend._
import cromwell.backend.callcaching.{BackendHashingMethods, DefaultBackendHashingMethods}
import wdl4s.Call
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}

class RetryableBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor, actorSystem: ActorSystem) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                serviceRegistryActor: ActorRef): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef): Props = {
    RetryableBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = NoFunctions

  override val backendHashingMethods: BackendHashingMethods = DefaultBackendHashingMethods(actorSystem)
}
