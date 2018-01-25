package cromwell.engine.workflow.lifecycle.execution.job.preparation

import cats.data.Validated.Valid
import cats.syntax.validated._
import cromwell.core.ExecutionIndex.ExecutionIndex
import common.validation.ErrorOr.ErrorOr
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import shapeless.Poly1
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OptionalInputDefinition}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

object InputPointerToWdlValue extends Poly1 {
  // Function that can transform any of the coproduct types to an ErrorOr[WomValue]
  type ToWdlValueFn = (Map[String, WomValue], IoFunctionSet, ValueStore, ExecutionIndex, InputDefinition) => ErrorOr[WomValue]

  implicit def fromWomValue: Case.Aux[WomValue, ToWdlValueFn] = at[WomValue] {
    womValue => (_: Map[String, WomValue], _: IoFunctionSet, _: ValueStore, _: ExecutionIndex, _: InputDefinition) =>
      womValue.validNel: ErrorOr[WomValue]
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWdlValueFn] = at[OutputPort] {
    port => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, valueStore: ValueStore, index: ExecutionIndex, inputDefinition: InputDefinition) =>
      (valueStore.resolve(index)(port), inputDefinition) match {
        case (v: Valid[WomValue], _) => v
        case (_, InputDefinitionWithDefault(_, _, defaultExpression)) =>
          evaluate(defaultExpression, knownValues, ioFunctions)
        case (_, OptionalInputDefinition(_, optionalType)) => optionalType.none.validNel
        case _ => s"Failed to lookup input value for required input ${port.name} at index $index in value store $valueStore".invalidNel
      }
  }

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWdlValueFn] = at[WomExpression] { 
    womExpression => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, _: ValueStore, _ : ExecutionIndex, _: InputDefinition) =>
      womExpression.evaluateValue(knownValues, ioFunctions): ErrorOr[WomValue]
  }

  def evaluate(womExpression: WomExpression, knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet): ErrorOr[WomValue] =
    womExpression.evaluateValue(knownValues, ioFunctions)
}
