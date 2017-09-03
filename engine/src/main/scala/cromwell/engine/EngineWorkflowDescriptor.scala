package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.CallCachingMode
import cromwell.core.path.PathBuilder
import wdl4s.wom.callable.WorkflowDefinition
import wdl4s.wom.graph.TaskCallNode

// TODO WOM: rename namespace to workflow
case class EngineWorkflowDescriptor(namespace: WorkflowDefinition,
                                    backendDescriptor: BackendWorkflowDescriptor,
                                    backendAssignments: Map[TaskCallNode, String],
                                    failureMode: WorkflowFailureMode,
                                    pathBuilders: List[PathBuilder],
                                    callCachingMode: CallCachingMode,
                                    parentWorkflow: Option[EngineWorkflowDescriptor] = None) {
  
  val rootWorkflow: EngineWorkflowDescriptor = parentWorkflow match {
    case Some(parent) => parent.rootWorkflow
    case None => this
  }

  def isRootWorkflow = rootWorkflow.parentWorkflow.isEmpty
  
  lazy val id = backendDescriptor.id
  lazy val workflow = backendDescriptor.workflow
  lazy val name = workflow.name
  lazy val knownValues = backendDescriptor.knownValues
  
  def getWorkflowOption(key: WorkflowOption) = backendDescriptor.getWorkflowOption(key)
}
