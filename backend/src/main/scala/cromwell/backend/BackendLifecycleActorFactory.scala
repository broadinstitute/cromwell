package cromwell.backend

import akka.actor.Props
import wdl4s.Call

trait BackendLifecycleActorFactory {
  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       calls: Seq[Call],
                                       configurationDescriptor: BackendConfigurationDescriptor): Option[Props]
  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor): Props
  def workflowFinalizationActorProps(): Option[Props]
}
