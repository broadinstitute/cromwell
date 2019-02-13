package wdl.transforms.base.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.{CommandPartElement, ExpressionElement, PlaceholderAttributeSet}
import wdl.model.draft3.graph.{ExpressionValueConsumer, GeneratedValueHandle}
import wdl.model.draft3.graph.expression._
import wdl.transforms.base.linking.graph.LinkedGraphMaker
import wdl.transforms.base.wdlom2wom.expression.WdlomWomExpression
import wom.callable.RuntimeEnvironment
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.LocalName
import wom.types.{WomArrayType, WomPrimitiveType, WomType}
import wom.values.{WomArray, WomBoolean, WomOptionalValue, WomPrimitive, WomValue}
import wom.{CommandPart, InstantiatedCommand}
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.transforms.base.linking.expression._
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl._

object CommandPartElementToWomCommandPart {
  def convert(commandPart: CommandPartElement,
              typeAliases: Map[String, WomType],
              availableHandles: Set[GeneratedValueHandle])
             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
              fileEvaluator: FileEvaluator[ExpressionElement],
              typeEvaluator: TypeEvaluator[ExpressionElement],
              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[CommandPart] = commandPart match {
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
  override def toString: String = stringCommandPartElement.value
  override def instantiate(inputsMap: Map[LocalName, WomValue],
                           functions: IoFunctionSet,
                           valueMapper: WomValue => WomValue,
                           runtimeEnvironment: RuntimeEnvironment): ErrorOr[List[InstantiatedCommand]] = List(InstantiatedCommand(stringCommandPartElement.value)).validNel
}

case class WdlomWomPlaceholderCommandPart(attributes: PlaceholderAttributeSet, expression: WdlomWomExpression) extends CommandPart {
  def attributesToString: String = attributes.toWdlV1
  // Yes, it's sad that we need to put ${} here, but otherwise we won't cache hit from draft-2 command sections
  override def toString: String = "${" + s"$attributesToString${expression.expressionElement.toWdlV1}" + "}"

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
      case None => InstantiatedCommand(commandString = "", createdFiles = value.sideEffectFiles.toList).validNel
    }
    case WomArray(WomArrayType(_ : WomPrimitiveType), arrayValue) => attributes.sepAttribute match {
      case Some(separator) => InstantiatedCommand(commandString = arrayValue.map(valueMapper(_).valueString).mkString(separator), createdFiles = value.sideEffectFiles.toList).validNel
      case None => "Array value was given but no 'sep' attribute was provided".invalidNel
    }
    case other => s"Cannot interpolate ${other.womType.stableName} into a command string with attribute set [$attributes]".invalidNel
  }
}
