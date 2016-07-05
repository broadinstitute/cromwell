package cromwell.engine.backend.mock

import akka.actor.Props
import cromwell.backend._
import wdl4s.Call
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}

class RetryableBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call]): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, initializationData: Option[BackendInitializationData]): Props = {
    RetryableBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = NoFunctions
}
