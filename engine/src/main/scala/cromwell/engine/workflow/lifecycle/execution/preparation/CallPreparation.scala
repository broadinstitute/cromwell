package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
import cats.data.NonEmptyList
import cats.data.Validated.Valid
import cromwell.backend.BackendJobDescriptor
import cromwell.core.CromwellGraphNode._
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import lenthall.validation.ErrorOr._
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.{WdlOptionalValue, WdlValue}
import wdl4s.wom.WomEvaluatedCallInputs
import wdl4s.wom.callable.Callable._
import wdl4s.wom.expression.{IoFunctionSet, WomExpression}

import scala.language.postfixOps

object CallPreparation {
  sealed trait CallPreparationActorCommands
  case class Start(outputStore: OutputStore) extends CallPreparationActorCommands

  trait CallPreparationActorResponse

  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends CallPreparationActorResponse

  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse

  // TODO WOM: This will be much simpler once wom gives us InputDefinition -> OutputPort / Expression mappings
  // Hopefully we shouldn't have to do any "string lookup" anymore
  def resolveAndEvaluateInputs(callKey: CallKey,
                               workflowDescriptor: EngineWorkflowDescriptor,
                               expressionLanguageFunctions: IoFunctionSet,
                               outputStore: OutputStore): ErrorOr[WomEvaluatedCallInputs] = {
    import cats.syntax.validated._
    import lenthall.validation.ErrorOr._
    import lenthall.validation.Validation._
    val call = callKey.scope

    val inputMappingsFromPreviousCalls: Map[String, WdlValue] = {
      call.inputPorts.flatMap({ inputPort =>
        val outputPort = inputPort.upstream
        // TODO WOM: scatters ?
        outputStore.get(outputPort, None) map { outputPort.fullyQualifiedName -> _ }
      }) toMap
    }

    // Previous calls outputs and known workflow level values
    val externalInputMappings = inputMappingsFromPreviousCalls ++ workflowDescriptor.backendDescriptor.knownValues

    def resolveInputDefinitionFold(accumulatedInputsSoFar: Map[InputDefinition, ErrorOr[WdlValue]],
                                   inputDefinition: InputDefinition): Map[InputDefinition, ErrorOr[WdlValue]] = {

      val validInputsAccumulated: Map[String, WdlValue] = accumulatedInputsSoFar.collect({
        case (input, Valid(errorOrWdlValue)) => input.name -> errorOrWdlValue
      }) ++ externalInputMappings

      // Tries to evaluate a WomExpression and coerce it to the desired type
      def evaluateAndCoerce(expression: WomExpression, coerceTo: WdlType): ErrorOr[WdlValue] = {
        expression.evaluateValue(validInputsAccumulated, expressionLanguageFunctions) flatMap {
          coerceTo.coerceRawValue(_).toErrorOr
        }
      }

      // Tries to find the current input definition in the map of known values
      def resolveFromKnownValues: ErrorOr[WdlValue] = {
        validInputsAccumulated.get(inputDefinition.name)
          // FullyQualified names are not working so make one up ourselves. Again this will be unnecessary when wom give us mappings
          .orElse(validInputsAccumulated.get(workflowDescriptor.workflow.name + "." + call.unqualifiedName + "." + inputDefinition.name)) match {
          case Some(value) => value.validNel
          case None => s"Can't find a known value for ${inputDefinition.name}".invalidNel
        }
      }

      // Tries to find the current input definition in the list of expression based inputs, and if so evaluates it
      def resolveAsExpression: ErrorOr[WdlValue] = {
        call.expressionBasedInputs.get(inputDefinition.name) map { instantiatedExpression =>
          evaluateAndCoerce(instantiatedExpression.expression, inputDefinition.womType)
        } getOrElse s"Can't find an expression to evaluate for ${inputDefinition.name}".invalidNel
      }

      // Tries to resolve the input definition by first looking in the known values, otherwise evaluate the associated
      // expression, if any. The order is important because we want to prioritize a value passed through workflow input over
      // an expression. This effectively allows to "override" the expression.
      def resolve: ErrorOr[WdlValue] = resolveFromKnownValues.orElse(resolveAsExpression)

      val evaluated = inputDefinition match {
        case OptionalInputDefinitionWithDefault(_, womType, default) => resolve.orElse(evaluateAndCoerce(default, womType))
        case OptionalInputDefinition(_, womType) => resolve.getOrElse(WdlOptionalValue.none(womType)).validNel
        case _: RequiredInputDefinition => resolve.leftMap(_ => NonEmptyList.of(s"Can't find an input value for ${inputDefinition.name}"))
      }

      accumulatedInputsSoFar + (inputDefinition -> evaluated)
    }

    call.callable.inputs.foldLeft(Map.empty[InputDefinition, ErrorOr[WdlValue]])(resolveInputDefinitionFold).sequence
  }
}
