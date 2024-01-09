package cromwell.engine

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.CallCachingMode
import cromwell.core.path.PathBuilder
import wom.callable.Callable
import wom.graph.CommandCallNode

case class EngineWorkflowDescriptor(topLevelCallable: Callable,
                                    backendDescriptor: BackendWorkflowDescriptor,
                                    backendAssignments: Map[CommandCallNode, String],
                                    failureMode: WorkflowFailureMode,
                                    pathBuilders: List[PathBuilder],
                                    callCachingMode: CallCachingMode,
                                    parentWorkflow: Option[EngineWorkflowDescriptor] = None
) {

  val rootWorkflow: EngineWorkflowDescriptor = parentWorkflow match {
    case Some(parent) => parent.rootWorkflow
    case None => this
  }

  def isRootWorkflow: Boolean = rootWorkflow.parentWorkflow.isEmpty

  lazy val id = backendDescriptor.id
  lazy val possiblyNotRootWorkflowId = id.toPossiblyNotRoot
  lazy val rootWorkflowId = rootWorkflow.id.toRoot
  lazy val callable = backendDescriptor.callable
  lazy val name = callable.name
  lazy val knownValues = backendDescriptor.knownValues

  def getWorkflowOption(key: WorkflowOption): Option[String] = backendDescriptor.getWorkflowOption(key)
}
