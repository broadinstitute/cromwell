package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, LoggingFSM}
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.Call

object WorkflowLifecycleActor {

  trait WorkflowLifecycleActorState { def terminal = false }
  trait WorkflowLifecycleActorTerminalState extends WorkflowLifecycleActorState { override val terminal = true }

  trait WorkflowLifecycleSuccessResponse extends EngineLifecycleStateCompleteResponse
  trait WorkflowLifecycleFailureResponse extends EngineLifecycleStateCompleteResponse

  object WorkflowLifecycleActorData {
    def empty = WorkflowLifecycleActorData(Map.empty, List.empty, Map.empty, List.empty)
  }

  /**
    * State data
    */
  case class WorkflowLifecycleActorData (backendActors: Map[ActorRef, String],
                                         successes: Seq[ActorRef],
                                         failures: Map[ActorRef, Throwable],
                                         aborted: Seq[ActorRef]) {

    def withBackendActors(backendActors: Map[ActorRef, String]) = this.copy(
      backendActors = this.backendActors ++ backendActors
    )
    def withSuccess(successfulActor: ActorRef) = this.copy(
      backendActors = this.backendActors - successfulActor,
      successes = successes :+ successfulActor)
    def withFailure(failedActor: ActorRef, reason: Throwable) = this.copy(
      backendActors = this.backendActors - failedActor,
      failures = failures + (failedActor -> reason))
    def withAborted(abortedActor: ActorRef) = this.copy(
      backendActors = this.backendActors - abortedActor,
      aborted = aborted :+ abortedActor
    )
  }
}

trait WorkflowLifecycleActor[S <: WorkflowLifecycleActorState] extends LoggingFSM[S, WorkflowLifecycleActorData] {

  val abortingState: S
  val successState: S
  val failureState: S
  val abortedState: S

  def successResponse: WorkflowLifecycleSuccessResponse
  def failureResponse(reasons: Seq[Throwable]): WorkflowLifecycleFailureResponse
  def abortedResponse: EngineLifecycleActorAbortedResponse

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> state if state.terminal =>
      log.info("State is now terminal. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      log.info(s"State is transitioning from $fromState to $toState.")
  }

  protected def checkForDoneAndTransition(newData: WorkflowLifecycleActorData): State = {
    if (checkForDone(newData)) {
      if (stateName == abortingState) {
        context.parent ! abortedResponse
        goto(abortedState) using newData
      } else {
        if (newData.failures.isEmpty) {
          context.parent ! successResponse
          goto(successState) using newData
        } else {
          context.parent ! failureResponse(newData.failures.values.toSeq)
          goto(failureState) using newData
        }
      }
    } else {
      stay using newData
    }
  }

  private def checkForDone(stateData: WorkflowLifecycleActorData) = stateData.backendActors.isEmpty
}
