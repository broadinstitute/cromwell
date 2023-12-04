package cromwell.engine.backend.mock

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import wom.expression.{IoFunctionSet, NoIoFunctionSet}
import wom.graph.CommandCallNode

import scala.concurrent.ExecutionContext

class RetryableBackendLifecycleActorFactory(val name: String,
                                            val configurationDescriptor: BackendConfigurationDescriptor
) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                ioActor: ActorRef,
                                                calls: Set[CommandCallNode],
                                                serviceRegistryActor: ActorRef,
                                                restarting: Boolean
  ): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      ioActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]
  ): Props =
    RetryableBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData],
                                           ioActorProxy: ActorRef,
                                           ec: ExecutionContext
  ): IoFunctionSet = NoIoFunctionSet
}
