package wom.expression

import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import common.validation.IOChecked.{IOChecked, _}
import shapeless.Poly1
import wom.callable.Callable.InputDefinition.InputValueMapper
import wom.callable.Callable.{InputDefinition, OverridableInputDefinitionWithDefault, OptionalInputDefinition}
import wom.graph.GraphNodePort.OutputPort
import wom.values.{WomOptionalValue, WomValue}

import scala.annotation.tailrec

object InputPointerToWomValue extends Poly1 {
  type OutputPortLookup = OutputPort => ErrorOr[WomValue]
  // Function that can transform any of the coproduct types to an ErrorOr[WomValue]
  type ToWomValueFn = (Map[String, WomValue], IoFunctionSet, OutputPortLookup, InputDefinition) => IOChecked[WomValue]

  def withInitializedValue(ioFunctionSet: IoFunctionSet, womValue: WomValue)(block: WomValue => IOChecked[WomValue]): IOChecked[WomValue] = {
    for {
      initialized <-  womValue.initialize(ioFunctionSet)
      evaluated <- block(initialized)
    } yield evaluated
  }

  implicit def fromWomValue: Case.Aux[WomValue, ToWomValueFn] = at[WomValue] {
    womValue => (_: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputPortLookup, inputDefinition: InputDefinition) =>
      withInitializedValue(ioFunctions, womValue) {
        inputDefinition.valueMapper(ioFunctions)
      }
  }

  implicit def fromOutputPort: Case.Aux[OutputPort, ToWomValueFn] = at[OutputPort] {
    port => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, outputPortLookup: OutputPortLookup, inputDefinition: InputDefinition) =>
      (outputPortLookup(port), inputDefinition) match {
        case (Valid(womValue), _) if isDefined(womValue) =>
          withInitializedValue(ioFunctions, womValue) {
            inputDefinition.valueMapper(ioFunctions)
          }
        case (_, OverridableInputDefinitionWithDefault(_, _, defaultExpression, valueMapper, _)) =>
          evaluateAndMap(defaultExpression, knownValues, valueMapper, ioFunctions)
        case (_, OptionalInputDefinition(_, optionalType, valueMapper, _)) => valueMapper(ioFunctions)(optionalType.none)
        case _ => s"Failed to lookup input value for required input ${port.internalName}".invalidIOChecked
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

  implicit def fromWomExpression: Case.Aux[WomExpression, ToWomValueFn] = at[WomExpression] {
    womExpression => (knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet, _: OutputPortLookup, inputDefinition: InputDefinition) =>
      evaluateAndMap(womExpression, knownValues, inputDefinition.valueMapper, ioFunctions)
  }

  def evaluate(womExpression: WomExpression, knownValues: Map[String, WomValue], ioFunctions: IoFunctionSet): ErrorOr[WomValue] =
    womExpression.evaluateValue(knownValues, ioFunctions)

  def evaluateAndMap(womExpression: WomExpression,
                     knownValues: Map[String, WomValue],
                     valueMapper: InputValueMapper,
                     ioFunctions: IoFunctionSet): IOChecked[WomValue] = for {
    evaluated <- evaluate(womExpression, knownValues, ioFunctions).toIOChecked
    mapped <- valueMapper(ioFunctions)(evaluated)
  } yield mapped
}
