package cromwell.engine

import java.nio.file.FileSystem

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingMode
import wdl4s._

final case class EngineWorkflowDescriptor(backendDescriptor: BackendWorkflowDescriptor,
                                          workflowInputs: WorkflowCoercedInputs,
                                          backendAssignments: Map[Call, String],
                                          failureMode: WorkflowFailureMode,
                                          engineFilesystems: List[FileSystem],
                                          callCachingMode: CallCachingMode) {
  def id = backendDescriptor.id
  def namespace = backendDescriptor.workflowNamespace
  def name = namespace.workflow.unqualifiedName
  def getWorkflowOption(key: WorkflowOption) = backendDescriptor.getWorkflowOption(key)
}
