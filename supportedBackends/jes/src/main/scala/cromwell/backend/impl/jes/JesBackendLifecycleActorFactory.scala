package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.core.{ExecutionStore, OutputStore}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

import scala.language.postfixOps

object JesBackendLifecycleActorFactory {
  implicit class Jessify(val genericInitializationData: Option[BackendInitializationData]) {
    def toJes = genericInitializationData collectFirst { case d: JesBackendInitializationData => d } get
  }
}

case class JesBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor) extends BackendLifecycleActorFactory {
  import JesBackendLifecycleActorFactory._

  val jesConfiguration = new JesConfiguration(configurationDescriptor)

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call]): Option[Props] = {
    Option(JesInitializationActor.props(workflowDescriptor, calls, jesConfiguration))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor, initializationData: Option[BackendInitializationData]): Props = {
    JesJobExecutionActor.props(jobDescriptor, jesConfiguration, initializationData.toJes)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call],
                                              executionStore: ExecutionStore,
                                              outputStore: OutputStore,
                                              initializationData: Option[BackendInitializationData]) = {
    Option(JesFinalizationActor.props(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore, initializationData.toJes))
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {

    val jesCallPaths = JesCallPaths(jobKey, workflowDescriptor, jesConfiguration, initializationData.toJes.backendFilesystem)
    new JesExpressionFunctions(List(initializationData.toJes.backendFilesystem), jesCallPaths.callContext)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    JesWorkflowPaths(workflowDescriptor, jesConfiguration, initializationData.toJes.backendFilesystem).rootPath
  }
}
