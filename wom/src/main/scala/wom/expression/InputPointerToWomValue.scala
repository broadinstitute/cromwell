package wom.expression

import cats.data.Validated.Valid
import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import shapeless.Poly1
import wom.callable.Callable.InputDefinition.InputValueMapper
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OptionalInputDefinition}
import wom.graph.GraphNodePort.OutputPort
import wom.values.{WomOptionalValue, WomValue}

import scala.annotation.tailrec

object InputPointerToWomValue extends Poly1 {
  type OutputPortLookup = OutputPort => ErrorOr[WomValue]
  // Function that can transform any of the coproduct types to an ErrorOr[WomValue]
  type ToWdlValueFn = (Map[String, WomValue], IoFunctionSet, OutputPortLookup, InputDefinition) => ErrorOr[WomValue]

  implicit def fromWomValue: Case.Aux[WomValue, ToWdlValueFn] = at[WomValue] {
    womValue => (_: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputPortLookup, inputDefinition: InputDefinition) =>
      inputDefinition.valueMapper(ioFunctions)(womValue)
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWdlValueFn] = at[OutputPort] {
    port => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, outputPortLookup: OutputPortLookup, inputDefinition: InputDefinition) =>
      (outputPortLookup(port), inputDefinition) match {
        case (Valid(womValue), _) if isDefined(womValue) =>
          inputDefinition.valueMapper(ioFunctions)(womValue)
        case (_, InputDefinitionWithDefault(_, _, defaultExpression, valueMapper)) =>
          evaluateAndMap(defaultExpression, knownValues, valueMapper, ioFunctions)
        case (_, OptionalInputDefinition(_, optionalType, valueMapper)) => valueMapper(ioFunctions)(optionalType.none)
        case _ => s"Failed to lookup input value for required input ${port.name}".invalidNel
      }
  }

  @tailrec
  private def isDefined(womValue: WomValue): Boolean = {
    womValue match {
      case WomOptionalValue(_, Some(innerWomValue)) => isDefined(innerWomValue)
      case WomOptionalValue(_, None) => false
      case _ => true
    }
  }

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWdlValueFn] = at[WomExpression] {
    womExpression => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputPortLookup, inputDefinition: InputDefinition) =>
      evaluateAndMap(womExpression, knownValues, inputDefinition.valueMapper, ioFunctions)
  }

  def evaluate(womExpression: WomExpression, knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet): ErrorOr[WomValue] =
    womExpression.evaluateValue(knownValues, ioFunctions)
  
  def evaluateAndMap(womExpression: WomExpression,
                     knownValues: Map[String, WomValue],
                     valueMapper: InputValueMapper,
                     ioFunctions: IoFunctionSet): ErrorOr[WomValue] = (for {
    evaluated <- evaluate(womExpression, knownValues, ioFunctions).toEither
    mapped <- valueMapper(ioFunctions)(evaluated).toEither
  } yield mapped).toValidated
}
