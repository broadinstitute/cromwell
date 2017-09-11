package cromwell.engine.backend.mock

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.core.NoIoFunctionSet
import wdl4s.wom.expression.IoFunctionSet
import wdl4s.wom.graph.TaskCallNode

class RetryableBackendLifecycleActorFactory(val name: String,
                                            val configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                ioActor: ActorRef,
                                                calls: Set[TaskCallNode],
                                                serviceRegistryActor: ActorRef,
                                                restarting: Boolean): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      ioActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    RetryableBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): IoFunctionSet = NoIoFunctionSet
}
