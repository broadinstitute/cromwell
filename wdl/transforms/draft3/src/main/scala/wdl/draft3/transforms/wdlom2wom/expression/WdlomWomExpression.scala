package wdl.draft3.transforms.wdlom2wom.expression

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.draft3.transforms.linking.expression.types._
import wdl.draft3.transforms.linking.expression.values._
import wdl.draft3.transforms.linking.expression.files._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types.WomType
import wom.values.{WomFile, WomValue}

final case class WdlomWomExpression(expressionElement: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]) extends WomExpression {
  override def sourceString: String = expressionElement.toString

  override def inputs: Set[String] = {
    expressionElement.expressionConsumedValueHooks map { hook => linkedValues(hook).linkableName }
  }

  def evaluateValueForPlaceholder(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: ForCommandInstantiationOptions): ErrorOr[EvaluatedValue[_]] =
    expressionElement.evaluateValue(inputValues, ioFunctionSet, Option(forCommandInstantiationOptions))

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
    expressionElement.evaluateValue(inputValues, ioFunctionSet, None) map { _.value }

  // NB types can be determined using the linked values, so we don't need the inputMap:
  override def evaluateType(inputMap: Map[String, WomType]): ErrorOr[WomType] = expressionElement.evaluateType(linkedValues)

  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
    expressionElement.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
}
