package cromwell.engine.workflow.lifecycle

import akka.actor.SupervisorStrategy.{Escalate, Stop}
import akka.actor.{ActorRef, LoggingFSM, _}
import cromwell.backend.BackendInitializationData
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

  case class BackendActorAndInitializationData(actor: ActorRef, data: Option[BackendInitializationData])

  /**
    * State data
    */
  case class WorkflowLifecycleActorData(actors: Set[ActorRef],
                                        successes: Seq[BackendActorAndInitializationData],
                                        failures: Map[ActorRef, Throwable],
                                        aborted: Seq[ActorRef]) {

    def withActors(actors: Set[ActorRef]) = this.copy(
      actors = this.actors ++ actors
    )
    def withSuccess(successfulActor: ActorRef, data: Option[BackendInitializationData] = None) = this.copy(
      actors = this.actors - successfulActor,
      successes = successes :+ BackendActorAndInitializationData(successfulActor, data))
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

  def successResponse(data: WorkflowLifecycleActorData): WorkflowLifecycleSuccessResponse
  def failureResponse(reasons: Seq[Throwable]): WorkflowLifecycleFailureResponse

  override def supervisorStrategy = AllForOneStrategy() {
    case ex: ActorInitializationException =>
      context.parent ! failureResponse(Seq(ex))
      context.stop(self)
      Stop
    case t => super.supervisorStrategy.decider.applyOrElse(t, (_: Any) => Escalate)
  }

  whenUnhandled {
    case unhandledMessage =>
      workflowLogger.warn(s"received an unhandled message: $unhandledMessage")
      stay
  }

  onTransition {
    case _ -> state if state.terminal =>
      workflowLogger.debug("State is now terminal. Stopping self.")
      context.stop(self)
    case fromState -> toState =>
      workflowLogger.debug(s"State is transitioning from $fromState to $toState.")
  }

  protected def checkForDoneAndTransition(newData: WorkflowLifecycleActorData): State = {
    if (checkForDone(newData)) {
      if (newData.failures.isEmpty) {
        context.parent ! successResponse(newData)
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
