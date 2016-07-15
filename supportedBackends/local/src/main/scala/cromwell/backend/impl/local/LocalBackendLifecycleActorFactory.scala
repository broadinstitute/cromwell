package cromwell.backend.impl.local

import java.util.concurrent.Executors

import akka.actor.Props
import cromwell.backend._
import cromwell.backend.io.{JobPaths, SharedFsExpressionFunctions}
import cromwell.core.CallContext
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions
import lenthall.config.ScalaConfig._

import scala.concurrent.ExecutionContext

case class LocalBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {

  private val poolSize = configurationDescriptor.backendConfig.getIntOr("pool-size", 10)
  private val ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(poolSize))

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call]): Option[Props] = {
    Option(LocalInitializationActor.props(workflowDescriptor, calls, configurationDescriptor).withDispatcher("akka.dispatchers.backend-dispatcher"))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, initializationData: Option[BackendInitializationData]): Props = {
    LocalJobExecutionActor.props(jobDescriptor, configurationDescriptor, ec).withDispatcher("akka.dispatchers.backend-dispatcher")
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey, None)
      val callContext = new CallContext(
        jobPaths.callRoot,
        jobPaths.stdout.toAbsolutePath.toString,
        jobPaths.stderr.toAbsolutePath.toString
      )

      new SharedFsExpressionFunctions(LocalJobExecutionActor.fileSystems, callContext)
  }
}
