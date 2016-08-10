package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Actor, ActorRef, Props}
import cromwell.backend._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{ExecutionStore, JobKey, OutputStore}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor._
import cromwell.webservice.WdlValueJsonFormatter
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.util.TryUtil
import wdl4s.values.WdlValue

import scala.util.{Failure, Success, Try}

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
            serviceRegistryActor: ActorRef) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new JobPreparationActor(executionData, jobKey, factory, initializationData, serviceRegistryActor))
  }
}

case class JobPreparationActor(executionData: WorkflowExecutionActorData,
                               jobKey: BackendJobDescriptorKey,
                               factory: BackendLifecycleActorFactory,
                               initializationData: Option[BackendInitializationData],
                               serviceRegistryActor: ActorRef)
  extends Actor with WdlLookup with WorkflowLogging {

  override val workflowDescriptor: EngineWorkflowDescriptor = executionData.workflowDescriptor
  override val workflowId = workflowDescriptor.id
  override val executionStore: ExecutionStore = executionData.executionStore
  override val outputStore: OutputStore = executionData.outputStore
  override val expressionLanguageFunctions = factory.expressionLanguageFunctions(
    workflowDescriptor.backendDescriptor, jobKey, initializationData)

  override def receive = {
    case Start =>
      val response = resolveAndEvaluateInputs(jobKey, expressionLanguageFunctions) flatMap { prepareJobExecutionActor } match {
        case Success(m) => m
        case Failure(f) => BackendJobPreparationFailed(jobKey, f)
      }

      context.parent ! response
      context stop self
    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }

  def prepareJobExecutionActor(inputs: Map[LocallyQualifiedName, WdlValue]): Try[BackendJobPreparationSucceeded] = {
    for {
      evaluatedRuntimeAttributes <- evaluateRuntimeAttributes(jobKey, expressionLanguageFunctions, inputs)
      jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, evaluatedRuntimeAttributes, inputs)
    } yield BackendJobPreparationSucceeded(jobDescriptor, factory.jobExecutionActorProps(jobDescriptor, initializationData, serviceRegistryActor))
  }

  // Split inputs map (= evaluated workflow declarations + coerced json inputs) into [init\.*].last
  private lazy val splitInputs = workflowDescriptor.backendDescriptor.inputs map { case (fqn, v) => splitFqn(fqn) -> v }

  def evaluateRuntimeAttributes(jobKey: BackendJobDescriptorKey,
                                wdlFunctions: WdlStandardLibraryFunctions,
                                evaluatedInputs: Map[LocallyQualifiedName, WdlValue]): Try[Map[String, WdlValue]] = {
    val tryInputs = evaluatedInputs map { case (x, y) => x -> Success(y) }
    val mapOfTries = jobKey.call.task.runtimeAttributes.attrs mapValues {
      expr => expr.evaluate(buildMapBasedLookup(tryInputs), wdlFunctions)
    }
    TryUtil.sequenceMap(mapOfTries).map(addDefaultsToAttributes(_))
  }

  def addDefaultsToAttributes(specifiedAttributes: Map[LocallyQualifiedName, WdlValue]): Map[LocallyQualifiedName, WdlValue] = {
    import WdlValueJsonFormatter._
    import spray.json._   // IGNORE INTELLIJ - this is deliberate...

    val workflowOptions = executionData.workflowDescriptor.backendDescriptor.workflowOptions
    def isUnspecifiedAttribute(name: String) = !specifiedAttributes.contains(name)

    val missing = factory.runtimeAttributeDefinitions filter { x => isUnspecifiedAttribute(x.name) }
    val defaults = missing map { x => (x, workflowOptions.getDefaultRuntimeOption(x.name)) } collect {
      case (runtimeAttributeDefinition, Success(jsValue)) => runtimeAttributeDefinition.name -> jsValue.convertTo[WdlValue]
      case (RuntimeAttributeDefinition(name, _, Some(defaultValue), _), _) => name -> defaultValue
    }

    specifiedAttributes ++ defaults
  }

  def resolveAndEvaluateInputs(jobKey: BackendJobDescriptorKey,
                               wdlFunctions: WdlStandardLibraryFunctions): Try[Map[LocallyQualifiedName, WdlValue]] = {
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
