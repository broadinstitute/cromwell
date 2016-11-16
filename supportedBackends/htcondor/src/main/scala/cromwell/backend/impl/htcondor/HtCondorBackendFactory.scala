package cromwell.backend.impl.htcondor

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend._
import cromwell.backend.impl.htcondor.caching.CacheActorFactory
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.core.{CallContext, WorkflowOptions}
import wdl4s.TaskCall
import wdl4s.expression.WdlStandardLibraryFunctions

import scala.util.{Failure, Success, Try}

case class HtCondorBackendFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory with StrictLogging {

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Set[TaskCall],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(HtCondorInitializationActor.props(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    HtCondorJobExecutionActor.props(jobDescriptor, configurationDescriptor, serviceRegistryActor, resolveCacheProviderProps(jobDescriptor.workflowDescriptor.workflowOptions))
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPathsWithDocker(jobKey, workflowDescriptor, configurationDescriptor.backendConfig)
    val callContext = CallContext(
      jobPaths.callExecutionRoot,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )

    new SharedFileSystemExpressionFunctions(HtCondorJobExecutionActor.pathBuilders, callContext)
  }

  private def resolveCacheProviderProps(workflowOptions: WorkflowOptions) = {
    val defaultCacheEnabled = configurationDescriptor.backendConfig.getBoolean("cache.enabled")
    val cacheEnabled: Boolean = getBooleanFromWfOptions(workflowOptions, "cacheEnabled", defaultCacheEnabled)

    if (cacheEnabled) {
      val defaultForceRewrite = configurationDescriptor.backendConfig.getBoolean("cache.forceRewrite")
      val cacheForceRewrite = getBooleanFromWfOptions(workflowOptions, "cacheForceRw", defaultForceRewrite)
      val provider = configurationDescriptor.backendConfig.getString("cache.provider")
      val cacheFactory = Class.forName(provider)
        .getConstructor(classOf[Config])
        .newInstance(configurationDescriptor.backendConfig)
        .asInstanceOf[CacheActorFactory]
      Option(cacheFactory.getCacheActorProps(cacheForceRewrite))
    }
    else None
  }

  private def getBooleanFromWfOptions(workflowOptions: WorkflowOptions, optionKey: String, defaultValue: Boolean) = {
    workflowOptions.get(optionKey) match {
      case Success(value) => Try(value.toBoolean).getOrElse {
        logger.warn(s"Could not get '$optionKey' attribute from workflow options. Falling back to default value.")
        defaultValue
      }
      case Failure(_) => defaultValue
    }
  }
}

