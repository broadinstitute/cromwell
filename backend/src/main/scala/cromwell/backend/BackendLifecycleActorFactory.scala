package cromwell.backend

import java.nio.file.Path

import akka.actor.Props
import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

trait BackendLifecycleActorFactory {
  def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                       calls: Seq[Call]): Option[Props]

  def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor): Props

  def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                     calls: Seq[Call],
                                     executionStore: BackendExecutionStore,
                                     outputStore: BackendOutputStore): Option[Props] = None

  def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                  jobKey: BackendJobDescriptorKey): WdlStandardLibraryFunctions

  def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config): Path = {
      new WorkflowPaths(workflowDescriptor, backendConfig).executionRoot
  }
}
