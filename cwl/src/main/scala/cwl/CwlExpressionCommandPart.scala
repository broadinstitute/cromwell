package cwl

import cats.data.NonEmptyList
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

import scala.util.Try
import cats.syntax.either._

import scala.language.postfixOps

case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] =

        ParameterContext().
          addLocalInputs(inputsMap).
          map(_.setRuntime(runtimeEnvironment)).
          flatMap(pc => expr.fold(EvaluateExpression).apply(pc).toEither.leftMap(e => NonEmptyList.one(e.getMessage))).
          map(_.valueString).
          map(s => InstantiatedCommand(s)).
          toValidated
}

// TODO: Dan to revisit making this an Either (and perhaps adding some other cases)
case class CommandLineBindingCommandPart(argument: CommandLineBinding) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: (WomValue) => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[InstantiatedCommand] =
    ParameterContext().
      addLocalInputs(inputsMap).
      flatMap{
        pc =>
          argument match {
            case CommandLineBinding(_, _, _, _, _, Some(StringOrExpression.Expression(expression)), Some(false)) =>
              expression.fold(EvaluateExpression).apply(pc).map(valueMapper).toEither.leftMap(e => NonEmptyList.one(e.getMessage))
            case CommandLineBinding(_, _, _, _, _, Some(StringOrExpression.String(string)), Some(false)) =>
              Right(WomString(string))
            // There's a fair few other cases to add, but until then...
            case other => throw new NotImplementedError(s"As-yet-unsupported command line binding: $other")
        }
      }.map(_.valueString).
      map(s => InstantiatedCommand(s)).
      toValidated
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
