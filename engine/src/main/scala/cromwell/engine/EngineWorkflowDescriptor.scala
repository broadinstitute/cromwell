package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import wdl4s._

final case class EngineWorkflowDescriptor(backendDescriptor: BackendWorkflowDescriptor,
                                          declarations: WorkflowCoercedInputs,
                                          backendAssignments: Map[Call, String],
                                          failureMode: WorkflowFailureMode) {
  def id = backendDescriptor.id
  def namespace = backendDescriptor.workflowNamespace
  def name = namespace.workflow.unqualifiedName
}
