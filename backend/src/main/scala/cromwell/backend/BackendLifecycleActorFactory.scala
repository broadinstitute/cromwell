package cromwell.backend

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend.callcaching.FileHashingActor
import cromwell.backend.callcaching.FileHashingActor.FileHashingFunction
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.wdl.OnlyPureFunctions
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.core.{ExecutionStore, OutputStore}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions


trait BackendLifecycleActorFactory {
  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       calls: Set[Call],
                                       serviceRegistryActor: ActorRef): Option[Props] = None

  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                             initializationData: Option[BackendInitializationData],
                             serviceRegistryActor: ActorRef,
                             backendSingletonActor: Option[ActorRef]): Props

  /**
    * Providing this method to generate Props for a cache hit copying actor is optional.
    * To implement it, add a function:
    * def cacheHitCopyingActorInner(jobDescriptor: BackendJobDescriptor,initializationData: Option[BackendInitializationData], serviceRegistryActor: ActorRef): Props
    * And then override this method to point to it:
    * override def cacheHitCopyingActorProps = Option(cacheHitCopyingActorInner _)
    *
    * Simples!
    */
  def cacheHitCopyingActorProps: Option[(BackendJobDescriptor, Option[BackendInitializationData], ActorRef) => Props] = None

  def backendSingletonActorProps: Option[Props] = None

  def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                     calls: Set[Call],
                                     executionStore: ExecutionStore,
                                     outputStore: OutputStore,
                                     initializationData: Option[BackendInitializationData]): Option[Props] = None

  def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                  jobKey: BackendJobDescriptorKey,
                                  initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = OnlyPureFunctions

  def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config, initializationData: Option[BackendInitializationData]): Path = {
    new WorkflowPaths(workflowDescriptor, backendConfig).executionRoot
  }

  def runtimeAttributeDefinitions(initializationDataOption: Option[BackendInitializationData]): Set[RuntimeAttributeDefinition] = Set.empty

  lazy val fileHashingFunction: Option[FileHashingFunction] = None
  lazy val fileHashingActorCount: Int = 50

  def fileHashingActorProps: Props = FileHashingActor.props(fileHashingFunction)

  def jobExecutionTokenType: JobExecutionTokenType = JobExecutionTokenType("Default", None)
}
