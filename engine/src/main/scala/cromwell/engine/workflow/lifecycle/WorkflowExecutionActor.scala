package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedResponse, BackendJobExecutionSucceededResponse, ExecuteJobCommand}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core.{WorkflowId, _}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.dummy.DummyBackendJobExecutionActor
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor._
import wdl4s.expression.WdlEvaluator.{WdlValueMapper, StringMapper}
import wdl4s.{LocallyQualifiedName, _}
import wdl4s.expression.{WdlEvaluator, WdlEvaluatorBuilder, WdlFunctions}
import wdl4s.values.WdlValue

import scala.language.postfixOps

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

  def startJob(call: Call) = {
    // TODO: Support indexes and retries:
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val evaluatorBuilder = new WdlEvaluatorBuilder(evaluatorFor(call))
    val inputs = resolveCallDeclarations(call)
    val jobDescriptor: BackendJobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, evaluatorBuilder, inputs)
    val executionActor = backendForExecution(jobDescriptor, workflowDescriptor.backendAssignments(call))
    executionActor ! ExecuteJobCommand
  }

  when(WorkflowExecutionPendingState) {
    case Event(StartExecutingWorkflowCommand, _) =>
      if (workflowDescriptor.namespace.workflow.calls.size == 1) {
        startJob(workflowDescriptor.namespace.workflow.calls.head)
        goto(WorkflowExecutionInProgressState)
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

  /**
    * Creates an appropriate BackendJobExecutionActor to run a single job, according to the backend assignments.
    */
  def backendForExecution(jobDescriptor: BackendJobDescriptor, backendName: String): ActorRef =
    if (backendName == "dummy") {
      // Don't judge me! This is obviously not "production ready"!
      val configDescriptor = BackendConfigurationDescriptor("", ConfigFactory.load())
      val props = DummyBackendJobExecutionActor.props(jobDescriptor, configDescriptor)
      val key = jobDescriptor.key
      context.actorOf(props, name = s"${jobDescriptor.descriptor.id}-BackendExecutionActor-${key.call.taskFqn}-${key.index}-${key.attempt}")
    } else {
      ??? //TODO: Implement!
    }


  // Split inputs map (= evaluated workflow declarations + coerced json inputs) into [init\.*].last
  private lazy val splitInputs = workflowDescriptor.backendDescriptor.inputs map {
    case (fqn, v) => fqn.splitFqn -> v
  }

  // Unqualified workflow level inputs
  private lazy val unqualifiedWorkflowInputs: Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == workflowDescriptor.namespace.workflow.unqualifiedName => inputName -> v
  }

  // Unqualified call inputs for a specific call, from the input json
  private def unqualifiedCallInputs(call: Call): Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == call.fullyQualifiedName => inputName -> v
  }

  // Workflow inputs + call inputs, potentially unevaluated
  private def staticInputsFor(call: Call): Map[LocallyQualifiedName, WdlValue] = {
    val resolvedDeclarations = resolveCallDeclarations(call) collect {
      case declaration => declaration.name -> declaration.unevaluatedValue
    }
    unqualifiedWorkflowInputs ++ (resolvedDeclarations toMap)
  }

  /**
    * Try to resolve call input declarations (= find the corresponding WdlExpression)
    * using workflow level inputs, input mappings, or declaration definition.
    * Note 1: workflow level declarations are evaluated WdlValues
    * Note 2: Unresolved declarations are removed from the list,
    * it's up to the backend to verify that all inputs are resolved if they want to.
    */
  private def resolveCallDeclarations(call: Call): Seq[ResolvedDeclaration] = {
    val inputsFromWorkflow = unqualifiedCallInputs(call)

    call.task.declarations collect {
      case decl if decl.expression.isDefined => decl.resolveWith(decl.expression.get)
      case decl if call.inputMappings.contains(decl.name) => decl.resolveWith(call.inputMappings(decl.name))
      case decl if inputsFromWorkflow.contains(decl.name) => decl.resolveWith(inputsFromWorkflow(decl.name))
    }
  }

  /**
    *
    * @param preValueMapper applies a String => String function to the identifier, before attempting to resolve it
    * @param postValueMapper applies a WdlValue => WdlValue function to the result of the evaluation
    * @return a function that takes engineFunctions and pre/post value mappers, and returns an Evaluator
    */
  private def evaluatorFor(call: Call)
                                 (engineFunctions: WdlFunctions[WdlValue],
                                  preValueMapper: StringMapper,
                                  postValueMapper: WdlValueMapper): WdlEvaluator = {
    // The lookup function will need to be aware of call outputs, shards, scatter variables etc...
    // This could be done by adding several type fo lookup (like in the old WA), or adding more information to the parameters map.
    def lookup = WdlExpression.standardLookupFunction(staticInputsFor(call), call.task.declarations, engineFunctions)

    new WdlEvaluator(lookup, engineFunctions, preValueMapper, postValueMapper)
  }
}
