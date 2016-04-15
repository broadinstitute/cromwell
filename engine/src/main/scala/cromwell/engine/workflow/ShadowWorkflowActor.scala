package cromwell.engine.workflow

import akka.actor.SupervisorStrategy.Escalate
import akka.actor._
import cromwell.core.{WorkflowOptions, WorkflowId}
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.WorkflowDescriptor
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor.{WorkflowExecutionFailedResponse, RestartExecutingWorkflowCommand, StartExecutingWorkflowCommand}
import cromwell.engine.workflow.lifecycle.WorkflowFinalizationActor.{WorkflowFinalizationSucceededResponse, WorkflowFinalizationFailedResponse, StartEngineFinalizationCommand}
import cromwell.engine.workflow.lifecycle.WorkflowInitializationActor.{WorkflowInitializationFailedResponse, WorkflowInitializationSucceededResponse, StartInitializationCommand}
import cromwell.engine.workflow.lifecycle.{WorkflowFinalizationActor, WorkflowExecutionActor, WorkflowInitializationActor}
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorFailure, MaterializeWorkflowDescriptorSuccess}
import cromwell.engine.workflow.ShadowWorkflowActor._
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowOutputs
import wdl4s.Call
import scala.util.Success
import scala.util.Failure

object ShadowWorkflowActor {

  /**
    * Commands to the ShadowWorkflowActor
    */
  sealed trait ShadowWorkflowActorCommand
  case object StartWorkflowCommand extends ShadowWorkflowActorCommand
  case object AbortWorkflowCommand extends ShadowWorkflowActorCommand

  /**
    * Responses from the ShadowWorkflowActor
    */
  sealed trait ShadowWorkflowActorResponse
  case class ShadowWorkflowSucceededResponse(outputs: WorkflowOutputs) extends ShadowWorkflowActorResponse
  case class ShadowWorkflowFailedResponse(reasons: Seq[Throwable]) extends Exception with ShadowWorkflowActorResponse

  /**
    * States for the ShadowWorkflowActor FSM
    */
  sealed trait ShadowWorkflowActorState { def terminal = false }
  sealed trait ShadowWorkflowActorTerminalState extends ShadowWorkflowActorState { override def terminal = true}

  /**
    * Waiting for a Start or Restart command.
    */
  case object WorkflowUnstartedState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor is created. It now needs to validate and create its WorkflowDescriptor
    */
  case object MaterializingWorkflowDescriptorState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor has a descriptor. It now needs to create its set of WorkflowBackendActors
    */
  case object InitializingWorkflowState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor has a descriptor and a set of backends. It now needs to execute the hell out of its jobs.
    */
  case object ExecutingWorkflowState extends ShadowWorkflowActorState

  /**
    * The WorkflowActor has completed. So we're now finalizing whatever needs to be finalized on the backends
    */
  case object FinalizingWorkflowState extends ShadowWorkflowActorState

  /**
    * We're done.
    */
  case object WorkflowSucceededState extends ShadowWorkflowActorTerminalState

  /**
    * We're done. But we were aborted in the middle of the lifecycle.
    */
  case object WorkflowAbortedState extends ShadowWorkflowActorTerminalState

  /**
    * We're done. But the workflow has failed.
    */
  case object WorkflowFailedState extends ShadowWorkflowActorTerminalState

  case class ShadowWorkflowActorData(backendAssignments: Map[Call, String])

  /**
    * Mode in which the workflow should be started:
    */
  sealed trait StartMode
  case object StartNewWorkflow extends StartMode
  case object RestartExistingWorkflow extends StartMode

  def props(workflowId: WorkflowId, startMode: StartMode, wdlSource: WorkflowSourceFiles): Props = Props(new ShadowWorkflowActor(workflowId, startMode, wdlSource))
}

/**
  * Class that orchestrates a single workflow... Shadow style.
  */
class ShadowWorkflowActor(workflowId: WorkflowId, startMode: StartMode, wdlSource: WorkflowSourceFiles) extends LoggingFSM[ShadowWorkflowActorState, ShadowWorkflowActorData] with ActorLogging{

  private val tag = s"$workflowId-${this.getClass.getSimpleName}"

  startWith(WorkflowUnstartedState, ShadowWorkflowActorData(Map.empty))

  override def supervisorStrategy: SupervisorStrategy = OneForOneStrategy() { case _ => Escalate }

  when(WorkflowUnstartedState) {
    case Event(StartWorkflowCommand, _) =>
      val materializeWorkflowDescriptorActor = context.actorOf(MaterializeWorkflowDescriptorActor.props(), name = s"MaterializeWorkflowDescriptorActor-$workflowId")
      materializeWorkflowDescriptorActor ! MaterializeWorkflowDescriptorActor.MaterializeWorkflow(workflowId, wdlSource)
      // Kill when done
      materializeWorkflowDescriptorActor ! PoisonPill
      goto(MaterializingWorkflowDescriptorState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(MaterializingWorkflowDescriptorState) {
    case Event(MaterializeWorkflowDescriptorSuccess(workflowDescriptor: WorkflowDescriptor), stateData) =>
      val newBackendAssignments: Map[Call, String] = assignCallsToBackends(workflowDescriptor.namespace.workflow.calls, workflowDescriptor.workflowOptions)
      val initializerActor = context.actorOf(WorkflowInitializationActor.props(tag, newBackendAssignments), name = s"EngineWorkflowInitializationActor-$workflowId")
      initializerActor ! StartInitializationCommand
      goto(InitializingWorkflowState) using stateData.copy(backendAssignments = newBackendAssignments)
    case Event(MaterializeWorkflowDescriptorFailure(reason: Throwable), _) =>
      // TODO: Handle workflow failure
      ???
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(InitializingWorkflowState) {
    case Event(WorkflowInitializationSucceededResponse, stateData) =>
      val executionActor = context.actorOf(WorkflowExecutionActor.props(tag, stateData.backendAssignments), name = s"EngineWorkflowExecutionActor-$workflowId")
      val commandToSend = startMode match {
        case StartNewWorkflow => StartExecutingWorkflowCommand
        case RestartExistingWorkflow => RestartExecutingWorkflowCommand
      }
      executionActor ! commandToSend
      goto(ExecutingWorkflowState)
    case Event(WorkflowInitializationFailedResponse(reasons), _) =>
      // TODO: Handle workflow failure
      ???
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(ExecutingWorkflowState) {
    case Event(WorkflowExecutionFailedResponse, stateData) =>
      val finalizationActor = context.actorOf(WorkflowFinalizationActor.props(tag, stateData.backendAssignments), name = s"EngineWorkflowFinalizationActor-$workflowId")
      finalizationActor ! StartEngineFinalizationCommand
      goto(FinalizingWorkflowState)
    case Event(WorkflowExecutionFailedResponse(reasons), _) =>
      // TODO: Handle workflow failure
      ???
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  when(FinalizingWorkflowState) {
    case Event(WorkflowFinalizationSucceededResponse, stateData) =>
      context.parent ! ShadowWorkflowSucceededResponse
      goto(WorkflowSucceededState)
    case Event(WorkflowFinalizationFailedResponse(reasons), _) =>
      context.parent ! ShadowWorkflowFailedResponse(reasons)
      goto(WorkflowFailedState)
    case Event(AbortWorkflowCommand, stateData) => ??? // TODO: Handle abort
  }

  onTransition {
    case _ -> newState if newState.terminal =>
      context.stop(self)
  }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message: $unhandledMessage")
      stay
  }

  private def assignCallsToBackends(calls: Seq[Call], workflowOptions: WorkflowOptions): Map[Call, String] = {
    val backendKey = "backend"
    // ------------------------------------------------------------------//
    // TODO: This should be obtained from the backend config factory. Part of issue #649
    // This is currently a stub. Should be removed as part of #649 (?)
    val listOfPluggedInBackends = Seq("local", "sge", "jes", "htcondor")
    val defaultBackend = "local"
    // ------------------------------------------------------------------//

    import scala.collection.breakOut
    val backendFromOptions = workflowOptions.get(backendKey)
    val callAssignmentMap = backendFromOptions match {
      case Success(backendName) =>
        // Use this backend for all the calls in the workflow
        calls.map(_ -> backendName)(breakOut): Map[Call, String]
      case Failure(reason) =>
        log.warning("Backend was not defined in the workflow options, trying to runtime-attributes based per task backend allocation.", reason)
        calls.map { call =>
          val runtimeAttributesMap = call.task.runtimeAttributes.attrs
          val backend = runtimeAttributesMap.get(backendKey) match {
            case Some(backendNameAsExp) =>
              // FIXME: [A Big one!] This expression should be interpreted as a String. How could we do that? Or get the literal value of the parser terminal?
              // This will always fallback to the default backend currently as it can't possibly hope to
              // find that expression in the Map of plugged in backends
              val backendNameAsString = backendNameAsExp.toWdlString
              if (listOfPluggedInBackends.find(_.equalsIgnoreCase(backendNameAsString)).isDefined) backendNameAsString else defaultBackend
            case None => defaultBackend
          }
          (call -> backend)
        }(breakOut): Map[Call, String]
    }
    log.debug(s"$tag: Call-to-Backend assignments: $callAssignmentMap")
    callAssignmentMap
  }
}
