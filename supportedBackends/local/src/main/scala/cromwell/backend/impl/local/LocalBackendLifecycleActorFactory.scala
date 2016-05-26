package cromwell.backend.impl.local

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.backend.io.{SharedFsExpressionFunctions, JobPaths}
import cromwell.core.CallContext
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

case class LocalBackendLifecycleActorFactory(config: Config) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                configurationDescriptor: BackendConfigurationDescriptor): Option[Props] = {
    Option(LocalInitializationActor.props(workflowDescriptor, calls, configurationDescriptor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      configurationDescriptor: BackendConfigurationDescriptor): Props = {
    LocalJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def workflowFinalizationActorProps(): Option[Props] = None

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           configurationDescriptor: BackendConfigurationDescriptor): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
      val callContext = new CallContext(
        jobPaths.callRoot,
        jobPaths.stdout.toAbsolutePath.toString,
        jobPaths.stderr.toAbsolutePath.toString
      )

      new SharedFsExpressionFunctions(LocalJobExecutionActor.fileSystems, callContext)
  }
}

