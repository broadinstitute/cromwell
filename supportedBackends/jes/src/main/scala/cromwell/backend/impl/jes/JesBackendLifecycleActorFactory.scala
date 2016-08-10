package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, ActorSystem, Props}
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.backend.impl.jes.callcaching.JesBackendHashingMethods
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.{ExecutionStore, OutputStore}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

import scala.language.postfixOps

object JesBackendLifecycleActorFactory {
  implicit class Jessify(val genericInitializationData: Option[BackendInitializationData]) {
    // This leaves the result in an `Option` as finalization will be called even if initialization has failed, and if
    // initialization fails there won't be any initialization data.  The various `.get`s that occur below are in instances
    // where the workflow has successfully gotten past initialization and the JES initialization data is defined.
    def toJes: Option[JesBackendInitializationData] = genericInitializationData collectFirst { case d: JesBackendInitializationData => d }
  }
}

case class JesBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor, actorSystem: ActorSystem) extends BackendLifecycleActorFactory {
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
    JesJobExecutionActor.props(jobDescriptor, jesConfiguration, initializationData.toJes.get, serviceRegistryActor).withDispatcher(BackendDispatcher)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call],
                                              executionStore: ExecutionStore,
                                              outputStore: OutputStore,
                                              initializationData: Option[BackendInitializationData]) = {
    Option(JesFinalizationActor.props(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore, initializationData.toJes).withDispatcher(BackendDispatcher))
  }

  override val backendHashingMethods = JesBackendHashingMethods(actorSystem)

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {

    val jesCallPaths = initializationData.toJes.get.workflowPaths.toJesCallPaths(jobKey)
    new JesExpressionFunctions(List(jesCallPaths.gcsFileSystemWithUserAuth), jesCallPaths.callContext)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData.toJes.get.workflowPaths.rootPath
  }
}
