package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.CallCachingMode
import cromwell.core.path.PathBuilder
import wdl4s._

final case class EngineWorkflowDescriptor(backendDescriptor: BackendWorkflowDescriptor,
                                          workflowInputs: WorkflowCoercedInputs,
                                          backendAssignments: Map[Call, String],
                                          failureMode: WorkflowFailureMode,
                                          pathBuilders: List[PathBuilder],
                                          callCachingMode: CallCachingMode) {
  def id = backendDescriptor.id
  def namespace = backendDescriptor.workflowNamespace
  def name = namespace.workflow.unqualifiedName
  def getWorkflowOption(key: WorkflowOption) = backendDescriptor.getWorkflowOption(key)
}
