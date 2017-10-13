package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
import cats.data.Validated.Valid
import cromwell.backend.BackendJobDescriptor
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import lenthall.validation.ErrorOr._
import lenthall.validation.Validation._
import wom.callable.Callable._
import wom.expression.IoFunctionSet
import wom.values.{WdlValue, WomEvaluatedCallInputs}

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
                               outputStore: OutputStore): ErrorOr[WomEvaluatedCallInputs] = {

    callKey.node.inputDefinitionMappings.foldLeft(Map.empty[InputDefinition, ErrorOr[WdlValue]]) {
      case (accumulatedInputsSoFar, (inputDefinition, pointer)) =>
        // We could have a lenthall method for this kind of "filtering valid values"
        val validInputsAccumulated: Map[String, WdlValue] = accumulatedInputsSoFar.collect({
          case (input, Valid(errorOrWdlValue)) => input.name -> errorOrWdlValue
        })
        
        val coercedValue = pointer.fold(InputPointerToWdlValue).apply(
          validInputsAccumulated, expressionLanguageFunctions, outputStore, callKey.index
        ) flatMap(inputDefinition.womType.coerceRawValue(_).toErrorOr)

        accumulatedInputsSoFar + (inputDefinition -> coercedValue)
    }.sequence
  }
}
