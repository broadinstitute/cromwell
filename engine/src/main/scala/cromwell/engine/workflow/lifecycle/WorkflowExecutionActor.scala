package cromwell.engine.workflow.lifecycle

import akka.actor.{FSM, LoggingFSM, Props}
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedResponse, BackendJobExecutionSucceededResponse, ExecuteJobCommand}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendLifecycleActorFactory}
import cromwell.core.WorkflowId
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackend}
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor._
import wdl4s._
import wdl4s.values.WdlValue

import scala.util.Success

object WorkflowExecutionActor {

  /**
    * States
    */
  sealed trait WorkflowExecutionActorState { def terminal = false }
  sealed trait WorkflowExecutionActorTerminalState extends WorkflowExecutionActorState { override val terminal = true }

  case object WorkflowExecutionPendingState extends WorkflowExecutionActorState
  case object WorkflowExecutionInProgressState extends WorkflowExecutionActorState
  case object WorkflowExecutionSuccessfulState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionFailedState extends WorkflowExecutionActorTerminalState
  case object WorkflowExecutionAbortedState extends WorkflowExecutionActorTerminalState

  /**
    * State data
    */
  final case class WorkflowExecutionActorData()

  /**
    * Commands
    */
  sealed trait WorkflowExecutionActorCommand
  case object StartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object RestartExecutingWorkflowCommand extends WorkflowExecutionActorCommand
  case object AbortExecutingWorkflowCommand extends WorkflowExecutionActorCommand

  /**
    * Responses
    */
  sealed trait WorkflowExecutionActorResponse
  case object WorkflowExecutionSucceededResponse extends WorkflowExecutionActorResponse
  case object WorkflowExecutionAbortedResponse extends WorkflowExecutionActorResponse
  final case class WorkflowExecutionFailedResponse(reasons: Seq[Throwable]) extends WorkflowExecutionActorResponse

  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor): Props = Props(WorkflowExecutionActor(workflowId, workflowDescriptor))
}

final case class WorkflowExecutionActor(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor) extends LoggingFSM[WorkflowExecutionActorState, WorkflowExecutionActorData] {

  val tag = self.path.name
  startWith(WorkflowExecutionPendingState, WorkflowExecutionActorData())

  /** PBE: the return value of WorkflowExecutionActorState is just temporary.
    *      This should probably return a Try[BackendJobDescriptor], Unit, Boolean,
    *      Try[ActorRef], or something to indicate if the job was started
    *      successfully.  Or, if it can fail to start, some indication of why it
    *      failed to start
    */
  private def startJob(jobKey: BackendJobDescriptorKey,
                       inputs: Map[FullyQualifiedName, WdlValue],
                       configDescriptor: BackendConfigurationDescriptor,
                       factory: BackendLifecycleActorFactory): WorkflowExecutionActorState = {
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, inputs)
    val jobExecutionActor = context.actorOf(
      factory.jobExecutionActorProps(
        jobDescriptor,
        BackendConfigurationDescriptor(configDescriptor.backendConfig, configDescriptor.globalConfig)
      ),
      s"$workflowId-BackendExecutionActor-${jobDescriptor.key.tag}"
    )
    jobExecutionActor ! ExecuteJobCommand
    WorkflowExecutionInProgressState
  }

  private def startJob(jobKey: BackendJobDescriptorKey, inputs: Map[FullyQualifiedName, WdlValue]): WorkflowExecutionActorState = {
    workflowDescriptor.backendAssignments.get(jobKey.call) match {
      case None =>
        val message = s"Could not start call ${jobKey.tag} because it was not assigned a backend"
        log.error(s"$tag $message")
        context.parent ! WorkflowExecutionFailedResponse(Seq(new Exception(message)))
        WorkflowExecutionFailedState
      case Some(backendName) =>
        val attemptedConfigurationDescriptor = BackendConfiguration.backendConfigurationDescriptor(backendName)
        val attemptedActorFactory = CromwellBackend.shadowBackendLifecycleFactory(backendName)

        (attemptedConfigurationDescriptor, attemptedActorFactory) match {
          case (Success(configDescriptor), Success(factory)) =>
            startJob(jobKey, inputs, configDescriptor, factory)
          case (_, _) =>
            val errors = List(
              attemptedActorFactory.failed.map(new Exception(s"Could not get BackendLifecycleActor for backend $backendName", _)).toOption,
              attemptedConfigurationDescriptor.failed.map(new Exception(s"Could not get BackendConfigurationDescriptor for backend $backendName", _)).toOption
            ).flatten
            context.parent ! WorkflowExecutionFailedResponse(errors)
            WorkflowExecutionFailedState
        }
    }
  }

  when(WorkflowExecutionPendingState) {
    case Event(StartExecutingWorkflowCommand, _) =>
      if (workflowDescriptor.namespace.workflow.calls.size == 1) {
        val jobKey = BackendJobDescriptorKey(workflowDescriptor.namespace.workflow.calls.head, None, 1)
        val nextState = startJob(jobKey, Map.empty)
        goto(nextState)
      } else {
        // TODO: We probably do want to support > 1 call in a workflow!
        sender ! WorkflowExecutionFailedResponse(Seq(new Exception("Execution is not implemented for call count != 1")))
        goto(WorkflowExecutionFailedState)
      }
    case Event(RestartExecutingWorkflowCommand, _) =>
      // TODO: Restart executing
      goto(WorkflowExecutionInProgressState)
    case Event(AbortExecutingWorkflowCommand, _) =>
      context.parent ! WorkflowExecutionAbortedResponse
      goto(WorkflowExecutionAbortedState)
  }

  when(WorkflowExecutionInProgressState) {
    case Event(BackendJobExecutionSucceededResponse(jobKey, callOutputs), stateData) =>
      log.info(s"Job ${jobKey.call.fullyQualifiedName} succeeded! Outputs: ${callOutputs.mkString("\n")}")
      context.parent ! WorkflowExecutionSucceededResponse
      goto(WorkflowExecutionSuccessfulState)
    case Event(BackendJobExecutionFailedResponse(jobKey, reason), stateData) =>
      log.warning(s"Job ${jobKey.call.fullyQualifiedName} failed! Reason: $reason")
      goto(WorkflowExecutionFailedState)
    case Event(AbortExecutingWorkflowCommand, stateData) => ??? // TODO: Implement!
    case Event(_, _) => ??? // TODO: Lots of extra stuff to include here...
  }

  when(WorkflowExecutionSuccessfulState) { FSM.NullFunction }
  when(WorkflowExecutionFailedState) { FSM.NullFunction }
  when(WorkflowExecutionAbortedState) { FSM.NullFunction }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"$tag received an unhandled message: $unhandledMessage in state: $stateName")
      stay
  }

  onTransition {
    case _ -> toState if toState.terminal =>
      log.info(s"$tag done. Shutting down.")
      context.stop(self)
    case fromState -> toState =>
      log.info(s"$tag transitioning from $fromState to $toState.")
  }
}
