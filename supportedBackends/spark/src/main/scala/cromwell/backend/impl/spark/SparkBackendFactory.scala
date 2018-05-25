package cromwell.backend.impl.spark

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.core.CallContext
import cromwell.core.path.DefaultPathBuilder
import wom.expression.IoFunctionSet
import wom.graph.CommandCallNode

import scala.concurrent.ExecutionContext

case class SparkBackendFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef,
                                                calls: Set[CommandCallNode], serviceRegistryActor: ActorRef, restarting: Boolean): Option[Props] = {
    Option(SparkInitializationActor.props(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      ioActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    SparkJobExecutionActor.props(jobDescriptor, configurationDescriptor, ioActor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData],
                                           ioActorProxy: ActorRef,
                                           ec: ExecutionContext): IoFunctionSet = {
    val jobPaths = JobPathsWithDocker(jobKey, workflowDescriptor, configurationDescriptor.backendConfig)
    val callContext = CallContext(
      jobPaths.callExecutionRoot,
      jobPaths.standardPaths,
      isDocker = true
    )

    new SharedFileSystemExpressionFunctions(SparkJobExecutionActor.DefaultPathBuilders, callContext, ioActorProxy, ec)
  }

  override def pathBuilders(initializationDataOption: Option[BackendInitializationData]) = List(DefaultPathBuilder)
}
