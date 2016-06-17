package cromwell.engine

import java.nio.file.FileSystem

import cromwell.backend.BackendWorkflowDescriptor
import wdl4s._

final case class EngineWorkflowDescriptor(backendDescriptor: BackendWorkflowDescriptor,
                                          workflowInputs: WorkflowCoercedInputs,
                                          backendAssignments: Map[Call, String],
                                          failureMode: WorkflowFailureMode,
                                          engineFilesystems: List[FileSystem]) {
  def id = backendDescriptor.id
  def namespace = backendDescriptor.workflowNamespace
  def name = namespace.workflow.unqualifiedName
  def getWorkflowOption(key: String) = backendDescriptor.workflowOptions.get(key).toOption
}
