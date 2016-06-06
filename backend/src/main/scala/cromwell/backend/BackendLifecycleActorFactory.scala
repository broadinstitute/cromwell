package cromwell.backend

import akka.actor.Props
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

trait BackendLifecycleActorFactory {
  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       calls: Seq[Call]): Option[Props]
  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor): Props
  def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                     calls: Seq[Call]): Option[Props]
  def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                  jobKey: BackendJobDescriptorKey): WdlStandardLibraryFunctions
}
