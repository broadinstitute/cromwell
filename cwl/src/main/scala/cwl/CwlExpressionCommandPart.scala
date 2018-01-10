package cwl

import cats.data.NonEmptyList
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cwl.CommandLineTool.CommandInputParameter
import cwl.command.ParentName
import mouse.all._
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values._
import wom.{CommandPart, InstantiatedCommand}

import scala.language.postfixOps
import scala.util.Try

case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {

    val pc = ParameterContext().
      addLocalInputs(inputsMap).
      setRuntime(runtimeEnvironment)

    val evaluatedExpression = expr.fold(EvaluateExpression).apply(pc).toEither.leftMap(e => NonEmptyList.one(e.getMessage))

    evaluatedExpression.map(_.valueString).map(InstantiatedCommand.apply(_)). toValidated
  }
}

// TODO: Dan to revisit making this an Either (and perhaps adding some other cases)
case class CommandLineBindingCommandPart(argument: CommandLineBinding) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] = {
    val pc = ParameterContext().addLocalInputs(inputsMap)

    val expressionEvaluationResult: Either[NonEmptyList[WorkflowStepInputId], WomValue] = (argument match {
      case CommandLineBinding(_, _, _, _, _, Some(StringOrExpression.Expression(expression)), Some(false)) =>
        expression.fold(EvaluateExpression).apply(pc).map(valueMapper).toEither.leftMap(e => NonEmptyList.one(e.getMessage))
      case CommandLineBinding(_, _, _, _, _, Some(StringOrExpression.String(string)), Some(false)) =>
        Right(WomString(string))
      // There's a fair few other cases to add, but until then...
      case other => throw new NotImplementedError(s"As-yet-unsupported command line binding: $other")
    })

    val commandString = expressionEvaluationResult.map(_.valueString)

    commandString.map(InstantiatedCommand.apply(_)).toValidated
  }
}

case class InputParameterCommandPart(commandInputParameter: CommandInputParameter)(implicit parentName: ParentName) extends CommandPart {

  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment) =
    Try {
      val womValue: WomValue = commandInputParameter match {
        case cip: CommandInputParameter =>
          val localizedId = LocalName(FullyQualifiedName(cip.id).id)
          inputsMap.get(localizedId).cata(
            valueMapper, throw new RuntimeException(s"could not find $localizedId in map $inputsMap"))

        // There's a fair few other cases to add, but until then...
        case other => throw new NotImplementedError(s"As-yet-unsupported commandPart from command input parameters: $other")
      }
      InstantiatedCommand(womValue.valueString)
    } toErrorOr
}
