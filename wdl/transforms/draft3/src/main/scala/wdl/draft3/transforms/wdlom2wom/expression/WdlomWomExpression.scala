package wdl.draft3.transforms.wdlom2wom.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.linking.expression.types._
import wdl.draft3.transforms.linking.expression.values._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types.WomType
import wom.values.{WomFile, WomValue}

final case class WdlomWomExpression(expressionElement: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]) extends WomExpression {
  override def sourceString: String = expressionElement.toString

  override def inputs: Set[String] = linkedValues.map(_._2.linkableName).toSet

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
    expressionElement.evaluateValue(inputValues, ioFunctionSet)

  // NB types can be determined using the linked values, so we don't need the inputMap:
  override def evaluateType(inputMap: Map[String, WomType]): ErrorOr[WomType] = expressionElement.evaluateType(linkedValues)

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {
    Set.empty[WomFile].validNel
  }
}
