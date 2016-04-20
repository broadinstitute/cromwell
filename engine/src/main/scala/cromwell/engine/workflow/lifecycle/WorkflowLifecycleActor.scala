package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, ActorRef, LoggingFSM}
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor._
import wdl4s.Call

object WorkflowLifecycleActor {

  trait WorkflowLifecycleActorState { def terminal = false }
  trait WorkflowLifecycleActorTerminalState extends WorkflowLifecycleActorState { override val terminal = true }

  trait WorkflowLifecycleSuccessResponse
  trait WorkflowLifecycleFailureResponse
  trait WorkflowLifecycleAbortedResponse

  /**
    * State data
    */
  trait WorkflowLifecycleActorData{
    /**
      * The mapping from ActorRef to the backend-specific actor for that backend
      */
    def backendActors: Map[ActorRef, String]

    /**
      * The list of successful backend-specific actors
      */
    def successes: Seq[ActorRef]

    /**
      * The list of failed backend-specific actors
      */
    def failures: Map[ActorRef, Throwable]
  }
}

trait WorkflowLifecycleActor[S <: WorkflowLifecycleActorState, D <: WorkflowLifecycleActorData] extends LoggingFSM[S, D] {

  val successState: S
  val failureState: S

  def successResponse: WorkflowLifecycleSuccessResponse
  def failureResponse(reasons: Seq[Throwable]): WorkflowLifecycleFailureResponse

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

  /**
    * Turns a mapping from Call => backendName into a mapping from backendName => List[Call]
    */
  protected def callAssignments(backendAssignments: Map[Call, String]): Map[String, List[Call]] =
    backendAssignments.groupBy(_._2).mapValues(_.keys.toList)

  /**
    * Based on backendAssignments, creates a BackendWorkflowInitializationActor per backend and maps the actor to the
    * name of the backend which it is initializing.
    */
  protected def backendWorkflowActors(backendAssignments: Map[Call, String]):  Map[ActorRef, String] = {
    val callAssignmentMap = callAssignments(backendAssignments)
    val backendsNeedingActors = backendAssignments.values.toSet
    backendsNeedingActors
      .map { backend => (backendActor(backend, callAssignmentMap(backend)), backend) } // Create the actors
      .collect { case (Some(actorRef), backend) => actorRef -> backend } // Only track the backends which have actors
      .toMap
  }

  /**
    * Makes an appropriate Backend actor for this backend. The call assignments are a list of calls which this
    * backend will be or has been requested to perform.
    */
  protected def backendActor(backendName: String, callAssignments: Seq[Call]): Option[ActorRef]

  protected def checkForDoneAndTransition(newData: D): State = {
    if (checkForDone(newData)) {
      if (newData.failures.isEmpty){
        context.parent ! successResponse
        goto(successState)
      }
      else {
        context.parent ! failureResponse(newData.failures.values.toSeq)
        goto(failureState) using newData
      }
    }
    else {
      stay using newData
    }
  }

  private def checkForDone(stateData: WorkflowLifecycleActorData) = stateData.backendActors.isEmpty
}
