package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, LoggingFSM}
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._

object WorkflowLifecycleActor {

  trait WorkflowLifecycleActorState { def terminal = false }
  trait WorkflowLifecycleActorTerminalState extends WorkflowLifecycleActorState { override val terminal = true }

  trait WorkflowLifecycleSuccessResponse extends EngineLifecycleStateCompleteResponse
  trait WorkflowLifecycleFailureResponse extends EngineLifecycleStateCompleteResponse

  object WorkflowLifecycleActorData {
    def empty = WorkflowLifecycleActorData(Set.empty, List.empty, Map.empty, List.empty)
  }

  /**
    * State data
    */
  case class WorkflowLifecycleActorData(actors: Set[ActorRef],
                                        successes: Seq[ActorRef],
                                        failures: Map[ActorRef, Throwable],
                                        aborted: Seq[ActorRef]) {

    def withActors(actors: Set[ActorRef]) = this.copy(
      actors = this.actors ++ actors
    )
    def withSuccess(successfulActor: ActorRef) = this.copy(
      actors = this.actors - successfulActor,
      successes = successes :+ successfulActor)
    def withFailure(failedActor: ActorRef, reason: Throwable) = this.copy(
      actors = this.actors - failedActor,
      failures = failures + (failedActor -> reason))
    def withAborted(abortedActor: ActorRef) = this.copy(
      actors = this.actors - abortedActor,
      aborted = aborted :+ abortedActor
    )
  }
}

trait AbortableWorkflowLifecycleActor[S <: WorkflowLifecycleActorState] extends WorkflowLifecycleActor[S] {
  val abortingState: S
  val abortedState: S

  def abortedResponse: EngineLifecycleActorAbortedResponse

  override protected def checkForDoneAndTransition(newData: WorkflowLifecycleActorData): State = {
    if (checkForDone(newData)) {
      if (stateName == abortingState) {
        context.parent ! abortedResponse
        goto(abortedState) using newData
      } else super.checkForDoneAndTransition(newData)
    } else {
      stay using newData
    }
  }
}

trait WorkflowLifecycleActor[S <: WorkflowLifecycleActorState] extends LoggingFSM[S, WorkflowLifecycleActorData] with WorkflowLogging {
  val successState: S
  val failureState: S

  def successResponse: WorkflowLifecycleSuccessResponse
  def failureResponse(reasons: Seq[Throwable]): WorkflowLifecycleFailureResponse


  whenUnhandled {
    case unhandledMessage =>
      workflowLogger.warn(s"received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> state if state.terminal =>
      workflowLogger.info("State is now terminal. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      workflowLogger.info(s"State is transitioning from $fromState to $toState.")
  }

  protected def checkForDoneAndTransition(newData: WorkflowLifecycleActorData): State = {
    if (checkForDone(newData)) {
      if (newData.failures.isEmpty) {
        context.parent ! successResponse
        goto(successState) using newData
      } else {
        context.parent ! failureResponse(newData.failures.values.toSeq)
        goto(failureState) using newData
      }
    } else {
      stay using newData
    }
  }

  protected final def checkForDone(stateData: WorkflowLifecycleActorData) = stateData.actors.isEmpty
}
