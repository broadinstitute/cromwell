package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CromwellGraphNode._
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.{OutputStore, WomOutputStore}
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.{WdlOptionalValue, WdlValue}
import wdl4s.wom.WomEvaluatedCallInputs
import wdl4s.wom.callable.Callable._
import wdl4s.wom.expression.{IoFunctionSet, WomExpression}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object CallPreparation {
  sealed trait CallPreparationActorCommands
  case class Start(outputStore: OutputStore) extends CallPreparationActorCommands

  trait CallPreparationActorResponse

  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends CallPreparationActorResponse

  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse

  def resolveAndEvaluateInputs(callKey: CallKey,
                               workflowDescriptor: EngineWorkflowDescriptor,
                               expressionLanguageFunctions: IoFunctionSet,
                               outputStore: OutputStore): Try[WomEvaluatedCallInputs] = {
    import cats.instances.list._
    import cats.syntax.traverse._
    import cats.syntax.validated._
    import lenthall.validation.ErrorOr._
    import lenthall.validation.Validation._
    val call = callKey.scope
    val womOutputStore: WomOutputStore = outputStore.toWomOutputStore

    val inputMappingsFromPreviousCalls: Map[String, WdlValue] = call.inputPorts collect {
      case inputPort if womOutputStore.get(inputPort.upstream).isDefined =>
        val outputPort = inputPort.upstream
        // TODO WOM: scatters ?
        // TODO WOM: clean up
        outputPort.fullyQualifiedName -> womOutputStore.get(outputPort).get
    } toMap

    // Previous calls outputs and known workflow level values
    val externalInputMappings = inputMappingsFromPreviousCalls ++ workflowDescriptor.backendDescriptor.knownValues

    def resolveInputDefinitionFold(accumulatedInputsSoFar: Map[InputDefinition, ErrorOr[WdlValue]],
                                   inputDefinition: InputDefinition): Map[InputDefinition, ErrorOr[WdlValue]] = {

      val validInputsAccumulated: Map[String, WdlValue] = accumulatedInputsSoFar.collect({
        case (input, Valid(errorOrWdlValue)) => input.name -> errorOrWdlValue
      }) ++ externalInputMappings

      def evaluateAndCoerce(expression: WomExpression, coerceTo: WdlType): ErrorOr[WdlValue] = {
        expression.evaluateValue(validInputsAccumulated, expressionLanguageFunctions) flatMap {
          coerceTo.coerceRawValue(_).toErrorOr
        }
      }
      
      def resolveFromKnownValues: ErrorOr[WdlValue] = validInputsAccumulated.get(inputDefinition.name) match {
        case Some(value) => value.validNel
        case None => s"Can't find a known value for ${inputDefinition.name}".invalidNel
      }

      def resolveAsExpression: ErrorOr[WdlValue] = {
        call.expressionBasedInputs.get(inputDefinition.name) map { instantiatedExpression =>
          evaluateAndCoerce(instantiatedExpression.expression, inputDefinition.womType)
        } getOrElse s"Can't find an expression to evaluate for ${inputDefinition.name}".invalidNel
      }

      def resolveAsInputPort: ErrorOr[WdlValue] = (for {
        inputPort <- call.portBasedInputs.find(_.name == inputDefinition.name)
        valueFromUpstream <- womOutputStore.get(inputPort.upstream)
      } yield valueFromUpstream) match {
        case Some(value) => value.validNel
        case None => s"Can't find upstream value for ${inputDefinition.name}".invalidNel
      }

      def resolve: ErrorOr[WdlValue] = resolveFromKnownValues.orElse(resolveAsInputPort).orElse(resolveAsExpression)

      val evaluated = inputDefinition match {
        case OptionalInputDefinitionWithDefault(_, womType, default) => resolve.orElse(evaluateAndCoerce(default, womType))
        case OptionalInputDefinition(_, womType) => resolve.getOrElse(WdlOptionalValue.none(womType)).validNel
        case _: RequiredInputDefinition => resolve.leftMap(_ => NonEmptyList.of(s"Can't find an input value for ${inputDefinition.name}"))
      }

      accumulatedInputsSoFar + (inputDefinition -> evaluated)
    }

    val evaluatedInputs: Map[InputDefinition, ErrorOr[WdlValue]] = {
      call.callable.inputs.foldLeft(Map.empty[InputDefinition, ErrorOr[WdlValue]])(resolveInputDefinitionFold)
    }

    // TODO WOM: cleanup
    evaluatedInputs.toList.traverse[ErrorOr, (InputDefinition, WdlValue)]({case (a, b) => b.map(c => (a,c))}).map(_.toMap) match {
      case Valid(res) => 
        Success(res)
      case Invalid(f) => Failure(new Exception(f.toList.mkString(", ")))
    }
  }
}
