package cromwell.backend.sfs

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory, BackendWorkflowDescriptor}
import cromwell.core.Dispatcher
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

import scala.concurrent.Promise

trait SharedFileSystemBackendLifecycleActorFactory extends BackendLifecycleActorFactory {

  def configurationDescriptor: BackendConfigurationDescriptor

  /**
    * Returns the main engine for async execution.
    *
    * @return the main engine for async execution.
    */
  def asyncJobExecutionActorClass: Class[_]

  /**
    * Returns true if the backend will generate special commands to run with docker.
    *
    * @return true if the backend will generate special commands to run with docker.
    */
  def supportsDocker = false

  def runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder = {
    SharedFileSystemValidatedRuntimeAttributesBuilder.default
  }

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call],
                                                serviceRegistryActor: ActorRef) = {
    val params = SharedFileSystemInitializationActorParams(serviceRegistryActor, workflowDescriptor,
      configurationDescriptor, calls, supportsDocker, runtimeAttributesBuilder)
    Option(Props(new SharedFileSystemInitializationActor(params)).withDispatcher(Dispatcher.BackendDispatcher))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationDataOption: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef) = {
    def propsCreator(completionPromise: Promise[BackendJobExecutionResponse]): Props = {
      val params = SharedFileSystemAsyncJobExecutionActorParams(serviceRegistryActor, jobDescriptor,
        configurationDescriptor, completionPromise, supportsDocker, runtimeAttributesBuilder, initializationDataOption)
      Props(asyncJobExecutionActorClass, params).withDispatcher(Dispatcher.BackendDispatcher)
    }

    Props(new SharedFileSystemJobExecutionActor(jobDescriptor, configurationDescriptor, propsCreator)).
      withDispatcher(Dispatcher.BackendDispatcher)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]):
  WdlStandardLibraryFunctions = {
    SharedFileSystemExpressionFunctions(workflowDescriptor, configurationDescriptor, jobKey, initializationData)
  }
}
