package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Actor, ActorRef, Props}
import cromwell.backend._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{ExecutionStore, JobKey, OutputStore}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor._
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.util.TryUtil
import wdl4s.values.WdlValue

import scala.util.{Failure, Success, Try}

final case class JobPreparationActor(executionData: WorkflowExecutionActorData,
                                     jobKey: BackendJobDescriptorKey,
                                     factory: BackendLifecycleActorFactory,
                                     initializationData: Option[BackendInitializationData],
                                     serviceRegistryActor: ActorRef,
                                     backendSingletonActor: Option[ActorRef])
  extends Actor with WdlLookup with WorkflowLogging {

  override lazy val workflowDescriptor: EngineWorkflowDescriptor = executionData.workflowDescriptor
  override lazy val workflowId = workflowDescriptor.id
  override lazy val executionStore: ExecutionStore = executionData.executionStore
  override lazy val outputStore: OutputStore = executionData.outputStore
  override lazy val expressionLanguageFunctions = factory.expressionLanguageFunctions(
    workflowDescriptor.backendDescriptor, jobKey, initializationData)

  override def receive = {
    case Start =>
      val response = resolveAndEvaluateInputs(jobKey, expressionLanguageFunctions) map { prepareJobExecutionActor }
      context.parent ! (response recover { case f => BackendJobPreparationFailed(jobKey, f) }).get
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }

  // Split inputs map (= evaluated workflow declarations + coerced json inputs) into [init\.*].last
  private lazy val splitInputs = workflowDescriptor.backendDescriptor.inputs map { case (fqn, v) => splitFqn(fqn) -> v }

  def resolveAndEvaluateInputs(jobKey: BackendJobDescriptorKey,
                               wdlFunctions: WdlStandardLibraryFunctions): Try[Map[LocallyQualifiedName, WdlValue]] = {
    import RuntimeAttributeDefinition.buildMapBasedLookup
    Try {
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

      TryUtil.sequenceMap(inputEvaluationAttempt, s"Input evaluation for Call ${call.fullyQualifiedName} failed")
    }.flatten
  }

  // Unqualified call inputs for a specific call, from the input json
  private def unqualifiedInputsFromInputFile(call: Call): Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == call.fullyQualifiedName => inputName -> v
  }

  private def prepareJobExecutionActor(inputEvaluation: Map[LocallyQualifiedName, WdlValue]): JobPreparationActorResponse = {
    import RuntimeAttributeDefinition.{addDefaultsToAttributes, evaluateRuntimeAttributes}
    val curriedAddDefaultsToAttributes = addDefaultsToAttributes(factory.runtimeAttributeDefinitions(initializationData), workflowDescriptor.backendDescriptor.workflowOptions) _

    (for {
      unevaluatedRuntimeAttributes <- Try(jobKey.call.task.runtimeAttributes)
      evaluatedRuntimeAttributes <- evaluateRuntimeAttributes(unevaluatedRuntimeAttributes, expressionLanguageFunctions, inputEvaluation)
      attributesWithDefault = curriedAddDefaultsToAttributes(evaluatedRuntimeAttributes)
      jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, attributesWithDefault, inputEvaluation)
    } yield BackendJobPreparationSucceeded(jobDescriptor, factory.jobExecutionActorProps(jobDescriptor, initializationData, serviceRegistryActor, backendSingletonActor))) match {
      case Success(s) => s
      case Failure(f) => BackendJobPreparationFailed(jobKey, f)
    }
  }
}

object JobPreparationActor {
  sealed trait JobPreparationActorCommands
  case object Start extends JobPreparationActorCommands

  sealed trait JobPreparationActorResponse
  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends JobPreparationActorResponse
  case class BackendJobPreparationFailed(jobKey: JobKey, throwable: Throwable) extends JobPreparationActorResponse

  def props(executionData: WorkflowExecutionActorData,
            jobKey: BackendJobDescriptorKey,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            serviceRegistryActor: ActorRef,
            backendSingletonActor: Option[ActorRef]) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new JobPreparationActor(executionData, jobKey, factory, initializationData, serviceRegistryActor, backendSingletonActor))
  }
}
