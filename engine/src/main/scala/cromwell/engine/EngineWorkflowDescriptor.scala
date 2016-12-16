package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.CallCachingMode
import cromwell.core.path.PathBuilder
import wdl4s._

final case class EngineWorkflowDescriptor(namespace: WdlNamespaceWithWorkflow,
                                          backendDescriptor: BackendWorkflowDescriptor,
                                          backendAssignments: Map[TaskCall, String],
                                          failureMode: WorkflowFailureMode,
                                          pathBuilders: List[PathBuilder],
                                          callCachingMode: CallCachingMode,
                                          parentWorkflow: Option[EngineWorkflowDescriptor] = None) {
  
  val rootWorkflow: EngineWorkflowDescriptor = parentWorkflow match {
    case Some(parent) => parent.rootWorkflow
    case None => this
  }
  
  lazy val id = backendDescriptor.id
  lazy val workflow = backendDescriptor.workflow
  lazy val name = workflow.unqualifiedName
  lazy val knownValues = backendDescriptor.knownValues
  
  def getWorkflowOption(key: WorkflowOption) = backendDescriptor.getWorkflowOption(key)
}
