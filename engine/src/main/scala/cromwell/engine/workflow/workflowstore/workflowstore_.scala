package cromwell.engine.workflow.workflowstore

import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreState.StartableState

sealed trait WorkflowStoreState {
  def stateName: String
  def restarted: Boolean = false

  /**
    * Return the restarted version of this state
    */
  def toRestartState: WorkflowStoreState
}

object WorkflowStoreState {
  private val RunningStateName = "Running"
  private val AbortingStateName = "Aborting"
  private val SubmittedStateName = "Submitted"
  
  val AllStates = List(Running, Aborting, Submitted, RestartableRunning, RestartableAborting)
  
  def fromNameAndRestarted(stateName: String, restarted: Boolean): Option[WorkflowStoreState] = {
    List(Running, Aborting, Submitted) collectFirst {
      case state if state.stateName.equalsIgnoreCase(stateName) =>
        if (restarted) state.toRestartState else state
    }
  }
  
  sealed trait StartableState extends WorkflowStoreState {
    /**
      * After being fetched, what state should the workflow be set to
      */
    def afterFetchedState: WorkflowStoreState
  }

  sealed trait RestartableState extends StartableState { override def restarted = true }

  case object Running extends WorkflowStoreState {
    override val stateName = RunningStateName
    override val toRestartState = RestartableRunning
  }

  case object Aborting extends WorkflowStoreState {
    override val stateName = AbortingStateName
    override val toRestartState = RestartableAborting
  }
  
  case object Submitted extends StartableState {
    override val stateName = SubmittedStateName
    override val toRestartState = this
    override val afterFetchedState = Running
  }

  /* 
    * Restartable states - Those states do not exist as is in the database,
    * they represent restartable workflows in Running and Aborting state. They exist so that we can type check that
    * a workflow is being started only if it's in a Startable state.
   */
  case object RestartableRunning extends RestartableState {
    override val stateName = RunningStateName
    override val toRestartState = this
    override val afterFetchedState = Running
  }

  case object RestartableAborting extends RestartableState {
    override val stateName = AbortingStateName
    override val toRestartState = this
    override val afterFetchedState = Aborting
  }
}

final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFilesCollection, state: StartableState)
