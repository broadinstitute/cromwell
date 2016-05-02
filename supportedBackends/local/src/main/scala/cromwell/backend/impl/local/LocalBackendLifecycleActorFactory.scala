package cromwell.backend.impl.local

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.{BackendWorkflowDescriptor, BackendConfigurationDescriptor, BackendJobDescriptor, BackendLifecycleActorFactory}
import wdl4s.Call

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
}

