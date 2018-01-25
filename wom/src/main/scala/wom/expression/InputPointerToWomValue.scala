package wom.expression

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OptionalInputDefinition}
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

object InputPointerToWomValue extends Poly1 {
  type OutputPortLookup = OutputPort => ErrorOr[WomValue]
  // Function that can transform any of the coproduct types to an ErrorOr[WomValue]
  type ToWdlValueFn = (Map[String, WomValue], IoFunctionSet, OutputPortLookup, InputDefinition) => ErrorOr[WomValue]

  implicit def fromWomValue: Case.Aux[WomValue, ToWdlValueFn] = at[WomValue] {
    womValue => (_: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputPortLookup, inputDefinition: InputDefinition) =>
      inputDefinition.valueMapper(ioFunctions)(womValue).validNel: ErrorOr[WomValue]
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWdlValueFn] = at[OutputPort] {
    port => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, outputPortLookup: OutputPortLookup, inputDefinition: InputDefinition) =>
      (outputPortLookup(port), inputDefinition) match {
        case (v: Valid[WomValue], _) => v.map(inputDefinition.valueMapper(ioFunctions))
        case (_, InputDefinitionWithDefault(_, _, defaultExpression, valueMapper)) =>
          evaluate(defaultExpression, knownValues, ioFunctions).map(valueMapper(ioFunctions))
        case (_, OptionalInputDefinition(_, optionalType, valueMapper)) => valueMapper(ioFunctions)(optionalType.none).validNel
        case _ => s"Failed to lookup input value for required input ${port.name}".invalidNel
      }
  }

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWdlValueFn] = at[WomExpression] {
    womExpression => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputPortLookup, inputDefinition: InputDefinition) =>
      womExpression.evaluateValue(knownValues, ioFunctions).map(inputDefinition.valueMapper(ioFunctions)): ErrorOr[WomValue]
  }

  def evaluate(womExpression: WomExpression, knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet): ErrorOr[WomValue] =
    womExpression.evaluateValue(knownValues, ioFunctions)
}
