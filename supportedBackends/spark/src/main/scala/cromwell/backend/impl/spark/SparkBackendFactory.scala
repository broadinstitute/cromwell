package cromwell.backend.impl.spark

import akka.actor.Props
import cromwell.backend._
import cromwell.backend.io.{SharedFsExpressionFunctions, JobPaths}
import cromwell.core.CallContext
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

case class SparkBackendFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call]): Option[Props] = {
    Option(SparkInitializationActor.props(workflowDescriptor, calls, configurationDescriptor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor): Props = {
    SparkJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
    val callContext = new CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )

    new SharedFsExpressionFunctions(SparkJobExecutionActor.fileSystems, callContext)
  }
}
