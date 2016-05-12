package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{FSM, LoggingFSM, Props}
import cromwell.backend._
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor._
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.util.TryUtil
import wdl4s.values.WdlValue

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

object JobPreparationActor {
  /** Commands */
  case object Start
  /** Responses */
  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, props: Props)
  case class BackendJobPreparationFailed(jobKey: JobKey, throwable: Throwable)
  /** States */
  sealed trait JobPreparationActorState
  case object Idle extends JobPreparationActorState
  case object EvaluatingInputs extends JobPreparationActorState
  case object BuildingProps extends JobPreparationActorState
  /** Internal Messages*/
  private case class EvaluatedInputs(inputs: Map[LocallyQualifiedName, WdlValue])
  private case class BackendActorProps(jobDescriptor: BackendJobDescriptor, props: Props)

  def props(executionData: WorkflowExecutionActorData,
            jobKey: BackendJobDescriptorKey,
            factory: BackendLifecycleActorFactory,
            configDescriptor: BackendConfigurationDescriptor) = Props(new JobPreparationActor(executionData, jobKey, factory, configDescriptor))
}

case class JobPreparationActor(executionData: WorkflowExecutionActorData,
                               jobKey: BackendJobDescriptorKey,
                               factory: BackendLifecycleActorFactory,
                               configDescriptor: BackendConfigurationDescriptor) extends LoggingFSM[JobPreparationActorState, Unit] with WdlLookup {

  override val workflowDescriptor: EngineWorkflowDescriptor = executionData.workflowDescriptor
  override val executionStore: ExecutionStore = executionData.executionStore
  override val outputStore: OutputStore = executionData.outputStore
  override val expressionLanguageFunctions = executionData.expressionLanguageFunctions

  startWith(Idle, Unit)

  when(Idle) {
    case Event(Start, _) =>
      resolveAndEvaluate(jobKey, expressionLanguageFunctions) onComplete {
        case Success(inputs) => self ! EvaluatedInputs(inputs)
        case Failure(f) =>
          context.parent ! BackendJobPreparationFailed(jobKey, f)
          context stop self
      }
      goto(EvaluatingInputs)
  }

  when(EvaluatingInputs) {
    case Event(EvaluatedInputs(inputs), _) =>
      val jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, inputs)
      val config = BackendConfigurationDescriptor(configDescriptor.backendConfig, configDescriptor.globalConfig)

      Future(factory.jobExecutionActorProps(jobDescriptor, config)) onComplete {
        case Success(props) =>
          context.parent ! BackendJobPreparationSucceeded(jobDescriptor, props)
          context stop self
        case Failure(f) =>
          context.parent ! BackendJobPreparationFailed(jobKey, f)
          context stop self
      }
      goto(BuildingProps)
  }

  when(BuildingProps) { FSM.NullFunction }

  whenUnhandled {
    case unhandledMessage =>
      log.warning(s"${self.path.name} received an unhandled message: $unhandledMessage in state: $stateName")
      stay
  }

  // Split inputs map (= evaluated workflow declarations + coerced json inputs) into [init\.*].last
  private lazy val splitInputs = workflowDescriptor.backendDescriptor.inputs map {
    case (fqn, v) => splitFqn(fqn) -> v
  }

  def resolveAndEvaluate(jobKey: BackendJobDescriptorKey,
                         wdlFunctions: WdlStandardLibraryFunctions): Future[Map[LocallyQualifiedName, WdlValue]] = Future {
    val call = jobKey.call
    lazy val callInputsFromFile = unqualifiedInputsFromInputFile(call)
    lazy val workflowScopedLookup = hierarchicalLookup(jobKey.call, jobKey.index) _

    // Try to resolve, evaluate and coerce declarations in order
    val inputEvaluationAttempt = call.task.declarations.foldLeft(Map.empty[LocallyQualifiedName, Try[WdlValue]])((inputs, declaration) => {
      val name = declaration.name

      // Try to resolve the declaration, and upon success evaluate the expression
      // If the declaration is resolved but can't be evaluated this will throw an evaluation exception
      // If it can't be resolved it's ignored and won't appear in the final input map
      val evaluated: Option[Try[WdlValue]] = declaration.expression match {
        // Static expression in the declaration
        case Some(expr) => Option(expr.evaluate(buildMapBasedLookup(inputs), wdlFunctions))
        // Expression found in the input mappings
        case None if call.inputMappings.contains(name) => Option(call.inputMappings(name).evaluate(workflowScopedLookup, wdlFunctions))
        // Expression found in the input file
        case None if callInputsFromFile.contains(name) => Option(Success(callInputsFromFile(name)))
        // Expression can't be found
        case _ => None
      }

      // Leave out unresolved declarations
      evaluated match {
        case Some(value) =>
          val coercedValue = value flatMap declaration.wdlType.coerceRawValue
          inputs + ((name, coercedValue))
        case None => inputs
      }
    })

    TryUtil.sequenceMap(inputEvaluationAttempt, s"Input evaluation for Call ${call.fullyQualifiedName} failed").get
  }

  // Unqualified call inputs for a specific call, from the input json
  private def unqualifiedInputsFromInputFile(call: Call): Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == call.fullyQualifiedName => inputName -> v
  }

  private def buildMapBasedLookup(evaluatedDeclarations: Map[LocallyQualifiedName, Try[WdlValue]])(identifier: String): WdlValue = {
    val successfulEvaluations = evaluatedDeclarations collect {
      case (k, v) if v.isSuccess => k -> v.get
    }
    successfulEvaluations.getOrElse(identifier, throw new WdlExpressionException(s"Could not resolve variable $identifier as a task input"))
  }

}
