package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.core.Dispatcher.BackendDispatcher
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
                                                calls: Seq[Call],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(JesInitializationActor.props(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor).withDispatcher(BackendDispatcher))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef): Props = {
    JesJobExecutionActor.props(jobDescriptor, jesConfiguration, initializationData.toJes, serviceRegistryActor).withDispatcher(BackendDispatcher)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call],
                                              executionStore: ExecutionStore,
                                              outputStore: OutputStore,
                                              initializationData: Option[BackendInitializationData]) = {
    Option(JesFinalizationActor.props(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore, initializationData.toJes).withDispatcher(BackendDispatcher))
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {

    val jesCallPaths = initializationData.toJes.workflowPaths.toJesCallPaths(jobKey)
    new JesExpressionFunctions(List(jesCallPaths.gcsFileSystemWithUserAuth), jesCallPaths.callContext)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData.toJes.workflowPaths.rootPath
  }
}
