package cromwell.backend.standard

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.backend.standard.callcaching._
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.path.Path
import cromwell.core.{CallOutputs, Dispatcher}
import wom.expression.IoFunctionSet
import wom.graph.CommandCallNode

import scala.concurrent.ExecutionContext

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
  lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[StandardInitializationActor]

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
  lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = Option(classOf[DefaultStandardCacheHitCopyingActor])

  /**
    * Returns the cache hit copying class.
    *
    * @return the cache hit copying class.
    */
  lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(classOf[DefaultStandardFileHashingActor])

  /**
    * Returns the finalization class.
    *
    * @return the finalization class.
    */
  lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] = Option(classOf[StandardFinalizationActor])

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                                serviceRegistryActor: ActorRef, restart: Boolean): Option[Props] = {
    val params = workflowInitializationActorParams(workflowDescriptor, ioActor, calls, serviceRegistryActor, restart)
    val props = Props(initializationActorClass, params).withDispatcher(Dispatcher.BackendDispatcher)
    Option(props)
  }

  def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                        serviceRegistryActor: ActorRef, restarting: Boolean): StandardInitializationActorParams = {
    DefaultInitializationActorParams(workflowDescriptor, ioActor, calls, serviceRegistryActor, configurationDescriptor, restarting)
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationDataOption: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      ioActor: ActorRef,
                                      backendSingletonActorOption: Option[ActorRef]): Props = {
    val params = jobExecutionActorParams(jobDescriptor, initializationDataOption, serviceRegistryActor,
      ioActor, backendSingletonActorOption)
    Props(new StandardSyncExecutionActor(params)).withDispatcher(Dispatcher.BackendDispatcher)
  }

  def jobExecutionActorParams(jobDescriptor: BackendJobDescriptor,
                              initializationDataOption: Option[BackendInitializationData],
                              serviceRegistryActor: ActorRef,
                              ioActor: ActorRef,
                              backendSingletonActorOption: Option[ActorRef]): StandardSyncExecutionActorParams = {
    DefaultStandardSyncExecutionActorParams(jobIdKey, serviceRegistryActor, ioActor, jobDescriptor, configurationDescriptor,
      initializationDataOption, backendSingletonActorOption, asyncExecutionActorClass, MinimumRuntimeSettings())
  }

  override def fileHashingActorProps:
  Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Option[ActorRef]) => Props] = {
    fileHashingActorClassOption map {
      standardFileHashingActor => fileHashingActorInner(standardFileHashingActor) _
    }
  }
  
  def fileHashingActorInner(standardFileHashingActor: Class[_ <: StandardFileHashingActor])
                               (jobDescriptor: BackendJobDescriptor,
                                initializationDataOption: Option[BackendInitializationData],
                                serviceRegistryActor: ActorRef,
                                ioActor: ActorRef,
                                fileHashCacheActor: Option[ActorRef]): Props = {
    val params = fileHashingActorParams(jobDescriptor, initializationDataOption, serviceRegistryActor, ioActor, fileHashCacheActor)
    Props(standardFileHashingActor, params).withDispatcher(BackendDispatcher)
  }

  def fileHashingActorParams(jobDescriptor: BackendJobDescriptor,
                             initializationDataOption: Option[BackendInitializationData],
                             serviceRegistryActor: ActorRef,
                             ioActor: ActorRef,
                             fileHashCacheActor: Option[ActorRef]): StandardFileHashingActorParams = {
    DefaultStandardFileHashingActorParams(
      jobDescriptor, initializationDataOption, serviceRegistryActor, ioActor, configurationDescriptor, fileHashCacheActor)
  }

  override def cacheHitCopyingActorProps:
  Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef, ActorRef, Int, Option[BlacklistCache]) => Props] = {
    cacheHitCopyingActorClassOption map {
      standardCacheHitCopyingActor => cacheHitCopyingActorInner(standardCacheHitCopyingActor) _
    }
  }

  def cacheHitCopyingActorInner(standardCacheHitCopyingActor: Class[_ <: StandardCacheHitCopyingActor])
                               (jobDescriptor: BackendJobDescriptor,
                                initializationDataOption: Option[BackendInitializationData],
                                serviceRegistryActor: ActorRef,
                                ioActor: ActorRef,
                                cacheCopyAttempt: Int,
                                blacklistCache: Option[BlacklistCache]): Props = {
    val params = cacheHitCopyingActorParams(jobDescriptor, initializationDataOption, serviceRegistryActor, ioActor, cacheCopyAttempt, blacklistCache)
    Props(standardCacheHitCopyingActor, params).withDispatcher(BackendDispatcher)
  }

  def cacheHitCopyingActorParams(jobDescriptor: BackendJobDescriptor,
                                 initializationDataOption: Option[BackendInitializationData],
                                 serviceRegistryActor: ActorRef,
                                 ioActor: ActorRef,
                                 cacheCopyAttempt: Int,
                                 blacklistCache: Option[BlacklistCache]): StandardCacheHitCopyingActorParams = {
    DefaultStandardCacheHitCopyingActorParams(
      jobDescriptor, initializationDataOption, serviceRegistryActor, ioActor, configurationDescriptor, cacheCopyAttempt, blacklistCache)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                              jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                              initializationData: Option[BackendInitializationData]): Option[Props] = {
    finalizationActorClassOption map { finalizationActorClass =>
      val params = workflowFinalizationActorParams(workflowDescriptor, ioActor, calls, jobExecutionMap, workflowOutputs,
        initializationData)
      Props(finalizationActorClass, params).withDispatcher(BackendDispatcher)
    }
  }

  def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                      jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                      initializationDataOption: Option[BackendInitializationData]):
  StandardFinalizationActorParams = {
    DefaultStandardFinalizationActorParams(workflowDescriptor, calls, jobExecutionMap, workflowOutputs,
      initializationDataOption, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationDataOption: Option[BackendInitializationData],
                                           ioActorProxy: ActorRef,
                                           ec: ExecutionContext):
  IoFunctionSet = {
    val standardInitializationData = BackendInitializationData.as[StandardInitializationData](initializationDataOption)
    val jobPaths = standardInitializationData.workflowPaths.toJobPaths(jobKey, workflowDescriptor)
    standardInitializationData.expressionFunctions(jobPaths, ioActorProxy, ec)
  }
  
  override def pathBuilders(initializationDataOption: Option[BackendInitializationData]) = {
    val standardInitializationData = BackendInitializationData.as[StandardInitializationData](initializationDataOption)
    standardInitializationData.workflowPaths.pathBuilders
  } 

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData match {
      case Some(data) => data.asInstanceOf[StandardInitializationData].workflowPaths.executionRoot
      case None => super.getExecutionRootPath(workflowDescriptor, backendConfig, initializationData)
    }
  }

  override def getWorkflowExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                            initializationData: Option[BackendInitializationData]): Path = {
    initializationData match {
      case Some(data) => data.asInstanceOf[StandardInitializationData].workflowPaths.workflowRoot
      case None => super.getWorkflowExecutionRootPath(workflowDescriptor, backendConfig, initializationData)
    }
  }

  override def runtimeAttributeDefinitions(initializationDataOption: Option[BackendInitializationData]): Set[RuntimeAttributeDefinition] = {
    val initializationData = BackendInitializationData.
      as[StandardInitializationData](initializationDataOption)

    initializationData.runtimeAttributesBuilder.definitions.toSet
  }
}
