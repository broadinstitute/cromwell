package wdl.transforms.base.wdlom2wom.expression

import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.expression._
import wdl.model.draft3.graph.{ExpressionValueConsumer, GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.expressionElementWriter
import wom.expression.{FileEvaluation, IoFunctionSet, WomExpression}
import wom.types._
import wom.values.WomValue

final case class WdlomWomExpression private (expressionElement: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                                            (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
                                             fileEvaluator: FileEvaluator[ExpressionElement],
                                             typeEvaluator: TypeEvaluator[ExpressionElement],
                                             valueEvaluator: ValueEvaluator[ExpressionElement]) extends WomExpression {
  override def sourceString: String = expressionElement.toWdlV1

  override def inputs: Set[String] = {
    expressionElement.expressionConsumedValueHooks map { hook => linkedValues(hook).linkableName }
  }

  def evaluateValueForPlaceholder(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: ForCommandInstantiationOptions): ErrorOr[EvaluatedValue[_]] =
    expressionElement.evaluateValue(inputValues, ioFunctionSet, Option(forCommandInstantiationOptions))

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
    expressionElement.evaluateValue(inputValues, ioFunctionSet, None) map { _.value }

  private lazy val evaluatedType = expressionElement.evaluateType(linkedValues)
  // NB types can be determined using the linked values, so we don't need the inputMap:
  override def evaluateType(inputMap: Map[String, WomType]): ErrorOr[WomType] = evaluatedType

  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = {
    expressionElement.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo) map { _ map {
      FileEvaluation(_, optional = areAllFileTypesInWomTypeOptional(coerceTo), secondary = false)
    }}
  }
}

object WdlomWomExpression {
  def make(expressionElement: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
          (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
           fileEvaluator: FileEvaluator[ExpressionElement],
           typeEvaluator: TypeEvaluator[ExpressionElement],
           valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[WdlomWomExpression] = {
    val candidate = WdlomWomExpression(expressionElement, linkedValues)
    candidate.evaluatedType.contextualizeErrors(s"process expression '${candidate.sourceString}'") map { _ => candidate }
  }
}
