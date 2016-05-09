package cromwell.backend

import akka.actor.Props
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

trait BackendLifecycleActorFactory {
  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       calls: Seq[Call],
                                       configurationDescriptor: BackendConfigurationDescriptor): Option[Props]
  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor): Props
  def workflowFinalizationActorProps(): Option[Props]
  def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                  jobKey: BackendJobDescriptorKey,
                                  configurationDescriptor: BackendConfigurationDescriptor): WdlStandardLibraryFunctions
}
