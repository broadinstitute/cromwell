package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.CommandPartElement
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.graph.{GeneratedValueHandle, LinkedGraph, UnlinkedConsumedValueHook}
import wom.callable.RuntimeEnvironment
import wom.expression.{IoFunctionSet, WomExpression}
import wom.{CommandPart, InstantiatedCommand}
import wom.graph.LocalName
import wom.values.WomValue
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.draft3.transforms.linking.expression._
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.draft3.transforms.linking.graph.LinkedGraphMaker
import wom.types.WomType

object CommandPartElementToWomCommandPart {
  def convert(commandPart: CommandPartElement, typeAliases: Map[String, WomType], availableHandles: Set[GeneratedValueHandle]): ErrorOr[CommandPart] = commandPart match {
    case s: StringCommandPartElement => WdlomWomStringCommandPart(s).validNel
    case p: PlaceholderCommandPartElement => {
      val attributes = Map.empty[String, String] // TODO: placeholder attributes
      val consumedValues = p.expressionElement.expressionConsumedValueHooks

      for {
        consumedValueLookup <- LinkedGraphMaker.makeConsumedValueLookup(consumedValues, availableHandles)
        womExpression <- p.expressionElement.makeWomExpression(typeAliases, consumedValueLookup)
                              _ <- validatePlaceholderType(womExpression, attributes)
      } yield WdlomWomPlaceholderCommandPart(attributes, womExpression)
    }
  }

  private def validatePlaceholderType(womExpression: WomExpression, attributes: Map[String, String]): ErrorOr[Unit] = ().validNel
}

case class WdlomWomStringCommandPart(stringCommandPartElement: StringCommandPartElement) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(stringCommandPartElement.value)).validNel
}

case class WdlomWomPlaceholderCommandPart(attributes: Map[String, String], expression: WomExpression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
    val inputsValues = inputsMap map { case (localName, value) => localName.value -> valueMapper(value) }
    expression.evaluateValue(inputsValues, functions) map { evaluatedExpression =>
      List(InstantiatedCommand(evaluatedExpression.valueString))
    }
  }
}
