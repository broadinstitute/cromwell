package cromwell.engine.backend.mock

import akka.actor.Props
import cromwell.backend._
import wdl4s.Call
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}

class RetryableBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {
  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call]): Option[Props] = None

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor): Props = {
    RetryableBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call]): Option[Props] = None

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey): WdlStandardLibraryFunctions = NoFunctions
}
