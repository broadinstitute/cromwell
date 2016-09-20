package cromwell.backend.impl.spark

import akka.actor.{ActorRef, ActorSystem, Props}
import cromwell.backend._
import cromwell.backend.sfs.SharedFileSystemExpressionFunctions
import cromwell.backend.io.JobPaths
import cromwell.core.CallContext
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

case class SparkBackendFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor, actorSystem: ActorSystem) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], serviceRegistryActor: ActorRef): Option[Props] = {
    Option(SparkInitializationActor.props(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef): Props = {
    SparkJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
    val callContext = new CallContext(
      jobPaths.callExecutionRoot,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )

    new SharedFileSystemExpressionFunctions(SparkJobExecutionActor.DefaultFileSystems, callContext)
  }
}
