package cromwell.backend.impl.jes

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory, BackendWorkflowDescriptor}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

case class JesBackendLifecycleActorFactory(config: Config) extends BackendLifecycleActorFactory {

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                configurationDescriptor: BackendConfigurationDescriptor): Option[Props] = {
    Option(JesInitializationActor.props(workflowDescriptor, calls, configurationDescriptor))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      configurationDescriptor: BackendConfigurationDescriptor): Props = {
    JesJobExecutionActor.props(jobDescriptor, configurationDescriptor)
  }

  override def workflowFinalizationActorProps(): Option[Props] = None

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           configurationDescriptor: BackendConfigurationDescriptor): WdlStandardLibraryFunctions = {

    val fileSystem = buildGcsFileSystem(configurationDescriptor, workflowDescriptor)
    val jesCallPaths = JesCallPaths(jobKey, fileSystem, workflowDescriptor, configurationDescriptor.backendConfig)
    new JesExpressionFunctions(List(fileSystem), jesCallPaths.callContext)
  }
}

