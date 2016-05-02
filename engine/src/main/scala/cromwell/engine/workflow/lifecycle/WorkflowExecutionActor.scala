package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionFailedResponse, BackendJobExecutionSucceededResponse, ExecuteJobCommand}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core.{EvaluatorBuilder, WorkflowId, _}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.dummy.DummyBackendJobExecutionActor
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor._
import wdl4s.expression.WdlFunctions
import wdl4s.values.WdlValue
import wdl4s.{Call, LocallyQualifiedName, WdlExpression}

import scala.language.postfixOps
import scala.util.{Success, Try}

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
    val evaluatorBuilder = new EvaluatorBuilder(evaluatorBuilderFor(call))
    val inputs = inputsFor(call)
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
    unqualifiedWorkflowInputs ++ inputsFor(call)
  }

  /**
    * Create a map of inputs for a specific call.
    * This map may contain WdlExpressions that will need to be evaluated by the backend.
    * TODO This can technically be done in the backend but seems very boilerplate and error prone. Maybe put it in a backend trait ?
    */
  private def inputsFor(call: Call): Map[LocallyQualifiedName, WdlValue] = {
    // Task declarations that have a static value assigned
    val staticDeclarations = call.task.declarations collect {
      case declaration if declaration.expression.isDefined => declaration.name -> declaration.expression.get
    } toMap

    staticDeclarations ++ unqualifiedCallInputs(call) ++ call.inputMappings
  }

  /**
    *
    * @param preValueMapper applies a String => String function to the identifier, before attempting to resolve it
    * @param postValueMapper applies a WdlValue => WdlValue function to the result of the evaluation
    * @return a function that takes engineFunctions and pre/post value mappers, and returns an Evaluator
    */
  private def evaluatorBuilderFor(call: Call)
                                 (engineFunctions: WdlFunctions[WdlValue],
                                  preValueMapper: StringMapper,
                                  postValueMapper: WdlValueMapper): Evaluator = {
    def lookup = {
      // The lookup function will need to be aware of call outputs, shards, scatter variables etc...
      // This could be done by adding several type fo lookup (like in the old WA), or adding more information to the parameters map.
      val standardLookupFunction = WdlExpression.standardLookupFunction(staticInputsFor(call), call.task.declarations, engineFunctions)
       standardLookupFunction compose preValueMapper
    }

    def evaluator(wdlValue: WdlValue) = {
      wdlValue match {
        case wdlExpression: WdlExpression => wdlExpression.evaluate(lookup, engineFunctions) map postValueMapper
        case v: WdlValue => Success(v)
      }
    }

    new Evaluator(evaluator)
  }
}
