package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.CallCachingMode
import cromwell.core.path.PathBuilder
import wdl4s._

final case class EngineWorkflowDescriptor(namespace: WdlNamespaceWithWorkflow,
                                          backendDescriptor: BackendWorkflowDescriptor,
                                          workflowInputs: WorkflowCoercedInputs,
                                          backendAssignments: Map[TaskCall, String],
                                          failureMode: WorkflowFailureMode,
                                          pathBuilders: List[PathBuilder],
                                          callCachingMode: CallCachingMode,
                                          parentWorkflow: Option[EngineWorkflowDescriptor] = None) {
  
  val rootWorkflow: EngineWorkflowDescriptor = parentWorkflow match {
    case Some(parent) => parent.rootWorkflow
    case None => this
  }
  
  val id = backendDescriptor.id
  lazy val workflow = backendDescriptor.workflow
  lazy val name = workflow.unqualifiedName
  val inputs = backendDescriptor.inputs
  def getWorkflowOption(key: WorkflowOption) = backendDescriptor.getWorkflowOption(key)
}
