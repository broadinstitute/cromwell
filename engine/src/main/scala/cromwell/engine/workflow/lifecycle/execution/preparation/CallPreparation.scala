package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
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
    
    val allInputMappings = inputMappingsFromPreviousCalls ++ workflowDescriptor.backendDescriptor.knownValues
    
    def resolveAsExpression(inputDefinition: InputDefinition): Option[ErrorOr[WdlValue]] = {
      call.expressionBasedInputs.get(inputDefinition.name) map { instantiatedExpression =>
        evaluateAndCoerce(instantiatedExpression.expression, inputDefinition.womType)
      }
    }
    
    def evaluateAndCoerce(expression: WomExpression, coerceTo: WdlType): ErrorOr[WdlValue] = {
      expression.evaluateValue(allInputMappings, expressionLanguageFunctions) flatMap {
        coerceTo.coerceRawValue(_).toErrorOr
      }
    }
    
    def resolveAsInputPort(inputDefinition: InputDefinition): Option[ErrorOr[WdlValue]] = {
      call.portBasedInputs.find(_.name == inputDefinition.name) map { inputPort =>
        womOutputStore.get(inputPort.upstream) match {
          case Some(value) => value.validNel
          case None => s"Can't find output value for ${inputPort.upstream}".invalidNel
        }
      }
    }
    
    def resolveFromInputs(inputDefinition: InputDefinition) = {
      resolveAsInputPort(inputDefinition).orElse(resolveAsExpression(inputDefinition))
    }
    
    def resolveInputDefinition(inputDefinition: InputDefinition): ErrorOr[WdlValue] = inputDefinition match {
      case DeclaredInputDefinition(_, womType, expression) => 
        evaluateAndCoerce(expression, womType)
      case optionalWithDefault @ OptionalInputDefinitionWithDefault(_, womType, default) =>
        resolveFromInputs(optionalWithDefault).getOrElse(evaluateAndCoerce(default, womType))
      case optional @ OptionalInputDefinition(_, womType) =>
        resolveFromInputs(optional).getOrElse(WdlOptionalValue.none(womType).validNel)
      case required: RequiredInputDefinition =>
        resolveFromInputs(required).getOrElse(s"Cannot find a value for ${inputDefinition.name}".invalidNel)
    }
    
    val evaluatedInputs: Map[InputDefinition, ErrorOr[WdlValue]] = (for {
      // TODO WOM: Should inputs be a Seq instead of Set ?
      // The issue is inputs can depend on each other, which means they need to be evaluated in "some" order
      // If it's a Set, order needs to be inferred and circular dependencies checked, a Seq would make it easier
      // Also WDL currently enforces the evaluation order based on the declaration order
      inputDefinition <- call.callable.inputs
    } yield inputDefinition -> resolveInputDefinition(inputDefinition)).toMap[InputDefinition, ErrorOr[WdlValue]]

    // TODO WOM: cleanup
    evaluatedInputs.toList.traverse[ErrorOr, (InputDefinition, WdlValue)]({case (a, b) => b.map(c => (a,c))}).map(_.toMap) match {
      case Valid(res) => Success(res)
      case Invalid(f) => Failure(new Exception(f.toList.mkString(", ")))
    }
  }
}
