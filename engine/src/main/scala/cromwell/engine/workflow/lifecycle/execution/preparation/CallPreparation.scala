package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
import cats.data.Validated.Valid
import cromwell.backend.BackendJobDescriptor
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.ValueStore
import common.validation.ErrorOr._
import common.validation.Validation._
import wom.callable.Callable._
import wom.expression.IoFunctionSet
import wom.values.{WomValue, WomEvaluatedCallInputs}

object CallPreparation {
  sealed trait CallPreparationActorCommands
  case class Start(valueStore: ValueStore) extends CallPreparationActorCommands

  trait CallPreparationActorResponse

  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends CallPreparationActorResponse

  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse

  def resolveAndEvaluateInputs(callKey: CallKey,
                               workflowDescriptor: EngineWorkflowDescriptor,
                               expressionLanguageFunctions: IoFunctionSet,
                               valueStore: ValueStore): ErrorOr[WomEvaluatedCallInputs] = {

    callKey.node.inputDefinitionMappings.foldLeft(Map.empty[InputDefinition, ErrorOr[WomValue]]) {
      case (accumulatedInputsSoFar, (inputDefinition, pointer)) =>
        // We could have a commons method for this kind of "filtering valid values"
        val validInputsAccumulated: Map[String, WomValue] = accumulatedInputsSoFar.collect({
          case (input, Valid(errorOrWdlValue)) => input.name -> errorOrWdlValue
        })
        
        val coercedValue = pointer.fold(InputPointerToWdlValue).apply(
          callKey.node, validInputsAccumulated, expressionLanguageFunctions, valueStore, callKey.index
        ) flatMap(inputDefinition.womType.coerceRawValue(_).toErrorOr)

        accumulatedInputsSoFar + (inputDefinition -> coercedValue)
    }.sequence
  }
}
