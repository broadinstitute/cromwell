package cromwell.engine.workflow.workflowstore

import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

sealed trait WorkflowStoreState {
  def isStartable: Boolean
  def isRestart: Boolean = false
}

object WorkflowStoreState {
  sealed trait StartableState extends WorkflowStoreState {
    override def isStartable = true 
    def afterFetchedState: WorkflowStoreState
  }
  sealed trait NonStartableState extends WorkflowStoreState { 
    override def isStartable = false 
  }
  sealed trait RestartedState extends WorkflowStoreState {
    override def isRestart = true
  }
  case object Running extends NonStartableState
  case object Aborting extends NonStartableState
  
  case object Submitted extends StartableState {
    def afterFetchedState: WorkflowStoreState = Running
  }
  
  case object RestartableRunning extends StartableState with RestartedState {
    def afterFetchedState: WorkflowStoreState = Running
  }
  case object RestartableAborting extends StartableState with RestartedState {
    def afterFetchedState: WorkflowStoreState = Aborting
  }
}

final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFilesCollection, state: StartableState)
