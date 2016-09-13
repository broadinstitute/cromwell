package cromwell.backend

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend.callcaching.FileHashingActor
import cromwell.backend.callcaching.FileHashingActor.FileHashingFunction
import cromwell.backend.io.WorkflowPaths
import cromwell.core.{ExecutionStore, OutputStore}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions


trait BackendLifecycleActorFactory {
  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       calls: Seq[Call],
                                       serviceRegistryActor: ActorRef): Option[Props]

  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                             initializationData: Option[BackendInitializationData],
                             serviceRegistryActor: ActorRef): Props

  def cacheHitCopyingActorProps(jobDescriptor: BackendJobDescriptor,
                                initializationData: Option[BackendInitializationData],
                                serviceRegistryActor: ActorRef): Option[Props] = None

  def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                     calls: Seq[Call],
                                     executionStore: ExecutionStore,
                                     outputStore: OutputStore,
                                     initializationData: Option[BackendInitializationData]): Option[Props] = None

  def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                  jobKey: BackendJobDescriptorKey,
                                  initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions

  def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config, initializationData: Option[BackendInitializationData]): Path = {
    new WorkflowPaths(workflowDescriptor, backendConfig).executionRoot
  }

  def runtimeAttributeDefinitions(initializationDataOption: Option[BackendInitializationData]): Set[RuntimeAttributeDefinition] = Set.empty

  lazy val fileHashingFunction: Option[FileHashingFunction] = None
  lazy val fileHashingActorCount: Int = 50

  def fileHashingActorProps: Props = FileHashingActor.props(fileHashingFunction)
}
