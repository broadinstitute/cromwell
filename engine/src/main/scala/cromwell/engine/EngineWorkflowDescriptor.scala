package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.CallCachingMode
import cromwell.core.path.PathBuilder
import wom.callable.WorkflowDefinition
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph.TaskCallNode

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
  /**
    * OutputPorts that could not be mapped to a WdlValue 
    */
  lazy val defaultExpressions: Map[OutputPort, WomExpression] = knownValues flatMap {
    case (port, resolvedInput) => resolvedInput.select[WomExpression] map { port -> _ }
  }
  
  def getWorkflowOption(key: WorkflowOption) = backendDescriptor.getWorkflowOption(key)
}
