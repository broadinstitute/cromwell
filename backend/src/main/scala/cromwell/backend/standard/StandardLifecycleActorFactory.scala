package cromwell.backend.standard

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.core.{CallOutputs, Dispatcher}
import cromwell.core.Dispatcher.BackendDispatcher
import wdl4s.TaskCall
import wdl4s.expression.WdlStandardLibraryFunctions

/**
  * May be extended for using the standard sync/async backend pattern.
  */
trait StandardLifecycleActorFactory extends BackendLifecycleActorFactory {
  /**
    * Config values for the backend, and a pointer to the global config.
    *
    * This is the single parameter passed into each factory during creation.
    *
    * @return The backend configuration.
    */
  def configurationDescriptor: BackendConfigurationDescriptor

  /**
    * Returns the key to use for storing and looking up the job id.
    *
    * @return the key to use for storing and looking up the job id.
    */
  def jobIdKey: String

  /**
    * Returns the asynchronous executor class.
    *
    * @return the asynchronous executor class.
    */
  def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor]

  /**
    * Returns the initialization class.
    *
    * @return the initialization class.
    */
  def initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[StandardInitializationActor]

  /**
    * Returns the synchronous executor class. By default using the standard sync executor should be sufficient for most
    * implementations.
    *
    * @return the synchronous executor class.
    */
  lazy val syncExecutionActorClass: Class[_ <: StandardSyncExecutionActor] = classOf[StandardSyncExecutionActor]

  /**
    * Returns the cache hit copying class.
    *
    * @return the cache hit copying class.
    */
  lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] =
    Option(classOf[StandardCacheHitCopyingActor])

  /**
    * Returns the finalization class.
    *
    * @return the finalization class.
    */
  lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[StandardFinalizationActor])

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    val params = workflowInitializationActorParams(workflowDescriptor, calls, serviceRegistryActor)
    val props = Props(initializationActorClass, params).withDispatcher(Dispatcher.BackendDispatcher)
    Option(props)
  }

  def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                        serviceRegistryActor: ActorRef): StandardInitializationActorParams = {
    DefaultInitializationActorParams(workflowDescriptor, calls, serviceRegistryActor, configurationDescriptor)
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationDataOption: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActorOption: Option[ActorRef]): Props = {
    val params = jobExecutionActorParams(jobDescriptor, initializationDataOption, serviceRegistryActor,
      backendSingletonActorOption)
    Props(new StandardSyncExecutionActor(params)).withDispatcher(Dispatcher.BackendDispatcher)
  }

  def jobExecutionActorParams(jobDescriptor: BackendJobDescriptor,
                              initializationDataOption: Option[BackendInitializationData],
                              serviceRegistryActor: ActorRef,
                              backendSingletonActorOption: Option[ActorRef]): StandardSyncExecutionActorParams = {
    DefaultStandardSyncExecutionActorParams(jobIdKey, serviceRegistryActor, jobDescriptor, configurationDescriptor,
      initializationDataOption, backendSingletonActorOption, asyncExecutionActorClass)
  }

  override def cacheHitCopyingActorProps:
  Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef) => Props] = {
    cacheHitCopyingActorClassOption map {
      standardCacheHitCopyingActor => cacheHitCopyingActorInner(standardCacheHitCopyingActor) _
    }
  }

  def cacheHitCopyingActorInner(standardCacheHitCopyingActor: Class[_ <: StandardCacheHitCopyingActor])
                               (jobDescriptor: BackendJobDescriptor,
                                initializationDataOption: Option[BackendInitializationData],
                                serviceRegistryActor: ActorRef): Props = {
    val params = cacheHitCopyingActorParams(jobDescriptor, initializationDataOption, serviceRegistryActor)
    Props(standardCacheHitCopyingActor, params).withDispatcher(BackendDispatcher)
  }

  def cacheHitCopyingActorParams(jobDescriptor: BackendJobDescriptor,
                                 initializationDataOption: Option[BackendInitializationData],
                                 serviceRegistryActor: ActorRef): StandardCacheHitCopyingActorParams = {
    DefaultStandardCacheHitCopyingActorParams(
      jobDescriptor, initializationDataOption, serviceRegistryActor, configurationDescriptor)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                              jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                              initializationData: Option[BackendInitializationData]): Option[Props] = {
    finalizationActorClassOption map { finalizationActorClass =>
      val params = workflowFinalizationActorParams(workflowDescriptor, calls, jobExecutionMap, workflowOutputs,
        initializationData)
      Props(finalizationActorClass, params).withDispatcher(BackendDispatcher)
    }
  }

  def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                      jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                      initializationDataOption: Option[BackendInitializationData]):
  StandardFinalizationActorParams = {
    DefaultStandardFinalizationActorParams(workflowDescriptor, calls, jobExecutionMap, workflowOutputs,
      initializationDataOption, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationDataOption: Option[BackendInitializationData]):
  WdlStandardLibraryFunctions = {
    val standardInitializationData = BackendInitializationData.as[StandardInitializationData](initializationDataOption)
    val jobPaths = standardInitializationData.workflowPaths.toJobPaths(jobKey, workflowDescriptor)
    standardInitializationData.expressionFunctions(jobPaths)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData.get.asInstanceOf[StandardInitializationData].workflowPaths.executionRoot
  }

  override def getWorkflowExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                            initializationData: Option[BackendInitializationData]): Path = {
    initializationData.get.asInstanceOf[StandardInitializationData].workflowPaths.workflowRoot
  }

  override def runtimeAttributeDefinitions(initializationDataOption: Option[BackendInitializationData]):
  Set[RuntimeAttributeDefinition] = {
    val initializationData = BackendInitializationData.
      as[StandardInitializationData](initializationDataOption)

    initializationData.runtimeAttributesBuilder.definitions.toSet
  }
}
