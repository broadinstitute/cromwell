package cwl

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.CommandLineTool.CommandInputParameter
import cwl.command.ParentName
import mouse.all._
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values._
import wom.{CommandPart, InstantiatedCommand}

case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {
    val stringKeyMap = inputsMap map {
      case (LocalName(localName), value) => localName -> value
    }
    val parameterContext = ParameterContext(
        inputs = stringKeyMap, runtimeOption = Option(runtimeEnvironment)
      )
    ExpressionEvaluator.eval(expr, parameterContext) map { womValue =>
      InstantiatedCommand(womValue.valueString)
    }
  }
}

// TODO: Dan to revisit making this an Either (and perhaps adding some other cases)
case class CommandLineBindingCommandPart(argument: CommandLineBinding) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {

    val stringInputsMap = inputsMap map {
      case (LocalName(localName), value) => localName -> valueMapper(value)
    }
    val parameterContext = ParameterContext(inputs = stringInputsMap, runtimeOption = Option(runtimeEnvironment))

    val womValueErrorOr: ErrorOr[WomValue] = argument match {
      case CommandLineBinding(_, _, _, _, _, Some(StringOrExpression.Expression(expression)), Some(false)) =>
        ExpressionEvaluator.eval(expression, parameterContext) map valueMapper
      case CommandLineBinding(_, _, _, _, _, Some(StringOrExpression.String(string)), Some(false)) =>
        WomString(string).valid
      // There's a fair few other cases to add, but until then...
      case other => throw new NotImplementedError(s"As-yet-unsupported command line binding: $other")
    }

    womValueErrorOr map { womValue =>
      InstantiatedCommand(womValue.valueString)
    }
  }
}

case class InputParameterCommandPart(commandInputParameter: CommandInputParameter)(implicit parentName: ParentName) extends CommandPart {

  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment) =
    validate {
      val womValue: WomValue = commandInputParameter match {
        case cip: CommandInputParameter =>
          val localizedId = LocalName(FullyQualifiedName(cip.id).id)
          inputsMap.get(localizedId).cata(
            valueMapper, throw new RuntimeException(s"could not find $localizedId in map $inputsMap"))

        // There's a fair few other cases to add, but until then...
        case other => throw new NotImplementedError(s"As-yet-unsupported commandPart from command input parameters: $other")
      }
      InstantiatedCommand(womValue.valueString)
    }
}
