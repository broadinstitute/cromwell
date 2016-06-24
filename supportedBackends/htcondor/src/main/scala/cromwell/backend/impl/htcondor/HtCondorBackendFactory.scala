package cromwell.backend.impl.htcondor

import akka.actor.Props
import cromwell.backend._
import cromwell.backend.impl.htcondor.caching.CacheActorFactory
import cromwell.backend.io.{JobPaths, SharedFsExpressionFunctions}
import cromwell.core.CallContext
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

case class HtCondorBackendFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call]): Option[Props] = {
    Option(HtCondorInitializationActor.props(workflowDescriptor, calls, configurationDescriptor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor): Props = {
    HtCondorJobExecutionActor.props(jobDescriptor, configurationDescriptor, resolveCacheProviderProps())
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
    val callContext = new CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )

    new SharedFsExpressionFunctions(HtCondorJobExecutionActor.fileSystems, callContext)
  }

  private def resolveCacheProviderProps() = {
    val cacheEnabled = configurationDescriptor.backendConfig.getBoolean("cache.enabled")

    if (cacheEnabled) {
      val provider = configurationDescriptor.backendConfig.getString("cache.provider")
      val cacheFactory = Class.forName(provider)
        .getConstructor(classOf[com.typesafe.config.Config])
        .newInstance(configurationDescriptor.backendConfig)
        .asInstanceOf[CacheActorFactory]
      Option(cacheFactory.getCacheActorProps())
    }
    else None
  }
}

