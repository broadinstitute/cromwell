package cromwell.backend.impl.htcondor


import akka.actor.Props
import cromwell.backend._
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
    HtCondorJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call]): Option[Props] = None

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
}

