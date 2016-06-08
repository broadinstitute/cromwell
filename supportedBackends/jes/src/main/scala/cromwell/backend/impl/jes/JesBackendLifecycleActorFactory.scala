package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.impl.jes.io._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory, BackendWorkflowDescriptor}
import cromwell.core.{ExecutionStore, OutputStore}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

case class JesBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {

  val jesConfiguration = new JesConfiguration(configurationDescriptor)

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call]): Option[Props] = {
    Option(JesInitializationActor.props(workflowDescriptor, calls, jesConfiguration))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor): Props = {
    JesJobExecutionActor.props(jobDescriptor, jesConfiguration)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call],
                                              executionStore: ExecutionStore,
                                              outputStore: OutputStore) = {
    Option(JesFinalizationActor.props(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore))
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey): WdlStandardLibraryFunctions = {

    val fileSystem = buildFilesystem(workflowDescriptor, jesConfiguration.jesAttributes.gcsFilesystemAuth, jesConfiguration.googleConfig)
    val jesCallPaths = JesCallPaths(jobKey, workflowDescriptor, jesConfiguration)
    new JesExpressionFunctions(List(fileSystem), jesCallPaths.callContext)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config): Path = {
    JesWorkflowPaths(workflowDescriptor, jesConfiguration).rootPath
  }
}
