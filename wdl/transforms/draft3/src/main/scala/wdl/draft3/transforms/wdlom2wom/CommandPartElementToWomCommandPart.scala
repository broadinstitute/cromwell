package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{CommandPartElement, PlaceholderAttributeSet}
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.graph.GeneratedValueHandle
import wom.callable.RuntimeEnvironment
import wom.expression.{IoFunctionSet, WomExpression}
import wom.{CommandPart, InstantiatedCommand}
import wom.graph.LocalName
import wom.values.{WomArray, WomBoolean, WomOptionalValue, WomPrimitive, WomValue}
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.draft3.transforms.linking.expression._
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.draft3.transforms.linking.graph.LinkedGraphMaker
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions}
import wom.types.{WomArrayType, WomPrimitiveType, WomType}

object CommandPartElementToWomCommandPart {
  def convert(commandPart: CommandPartElement, typeAliases: Map[String, WomType], availableHandles: Set[GeneratedValueHandle]): ErrorOr[CommandPart] = commandPart match {
    case s: StringCommandPartElement => WdlomWomStringCommandPart(s).validNel
    case p: PlaceholderCommandPartElement => {
      val attributes = p.attributes
      val consumedValues = p.expressionElement.expressionConsumedValueHooks

      for {
        consumedValueLookup <- LinkedGraphMaker.makeConsumedValueLookup(consumedValues, availableHandles)
        womExpression <- p.expressionElement.makeWomExpression(typeAliases, consumedValueLookup)
        _ <- validatePlaceholderType(womExpression, attributes)
      } yield WdlomWomPlaceholderCommandPart(attributes, womExpression.asInstanceOf[WdlomWomExpression])
    }
  }

  private def validatePlaceholderType(womExpression: WomExpression, attributes: PlaceholderAttributeSet): ErrorOr[Unit] = ().validNel
}

case class WdlomWomStringCommandPart(stringCommandPartElement: StringCommandPartElement) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(stringCommandPartElement.value)).validNel
}

case class WdlomWomPlaceholderCommandPart(attributes: PlaceholderAttributeSet, expression: WdlomWomExpression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = {
    val inputsValues = inputsMap map { case (localName, value) => localName.value -> value }
    expression.evaluateValueForPlaceholder(inputsValues, functions, ForCommandInstantiationOptions(valueMapper)) flatMap { evaluatedExpression =>
      instantiateFromValue(evaluatedExpression, valueMapper).map(List(_))
    }
  }

  private def instantiateFromValue(value: EvaluatedValue[_], valueMapper: WomValue => WomValue): ErrorOr[InstantiatedCommand] = value.value match {
    case WomBoolean(b) if attributes.trueAttribute.isDefined && attributes.falseAttribute.isDefined =>
      InstantiatedCommand(
        commandString = if (b) { attributes.trueAttribute.get } else { attributes.falseAttribute.get },
        createdFiles = value.sideEffectFiles.toList
      ).validNel
    case p: WomPrimitive => InstantiatedCommand(valueMapper(p).valueString, createdFiles = value.sideEffectFiles.toList).validNel
    case WomOptionalValue(_, Some(v)) => instantiateFromValue(value.copy(value = v), valueMapper)
    case WomOptionalValue(_, None) => attributes.defaultAttribute match {
      case Some(default) => InstantiatedCommand(commandString = default, createdFiles = value.sideEffectFiles.toList).validNel
      case None => "Optional value was not set and no 'default' attribute was provided".invalidNel
    }
    case WomArray(WomArrayType(_ : WomPrimitiveType), arrayValue) => attributes.sepAttribute match {
      case Some(separator) => InstantiatedCommand(commandString = arrayValue.map(valueMapper(_).valueString).mkString(separator), createdFiles = value.sideEffectFiles.toList).validNel
      case None => "Array value was given but no 'sep' attribute was provided".invalidNel
    }
    case other => s"Cannot interpolate ${other.womType.toDisplayString} into a command string with attribute set [$attributes]".invalidNel
  }
}
