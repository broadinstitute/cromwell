package cromwell.backend.impl.local

import java.util.concurrent.Executors

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.backend.io.{JobPaths, SharedFsExpressionFunctions}
import cromwell.core.CallContext
import cromwell.core.Dispatcher.BackendDispatcher
import lenthall.config.ScalaConfig._
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

import scala.concurrent.ExecutionContext
import scala.language.postfixOps

object LocalBackendLifecycleActorFactory {
  implicit class Localify(val genericInitializationData: Option[BackendInitializationData]) {
    def toLocal = genericInitializationData collectFirst { case d: LocalBackendInitializationData => d } get
  }
}

case class LocalBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {
  import LocalBackendLifecycleActorFactory._

  private val poolSize = configurationDescriptor.backendConfig.getIntOr("pool-size", 10)
  private val ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(poolSize))
  private val localConfiguration = new LocalConfiguration(configurationDescriptor)

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(LocalInitializationActor.props(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor, localConfiguration).withDispatcher(BackendDispatcher))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef): Props = {
    LocalJobExecutionActor.props(jobDescriptor, configurationDescriptor, serviceRegistryActor, initializationData.toLocal, ec).withDispatcher(BackendDispatcher)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
      val callContext = CallContext(
        jobPaths.callRoot,
        jobPaths.stdout.toAbsolutePath.toString,
        jobPaths.stderr.toAbsolutePath.toString
      )

      new SharedFsExpressionFunctions(initializationData.toLocal.workflowPaths.fileSystems, callContext)
  }
}
