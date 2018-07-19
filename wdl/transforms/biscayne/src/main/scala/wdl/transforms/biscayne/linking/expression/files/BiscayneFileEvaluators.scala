package wdl.transforms.biscayne.linking.expression.files

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{AsMap, AsPairs, CollectByKey}
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}
import wdl.transforms.biscayne.linking.expression.values.expressionEvaluator

object BiscayneFileEvaluators {
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] = new FileEvaluator[AsMap] {
    override def predictFilesNeededToEvaluate(a: AsMap, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {

      fileEvaluator.evaluateFilesNeededToEvaluate(a.param, inputs, ioFunctionSet, coerceTo)(fileEvaluator, expressionEvaluator)
    }
  }

  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] = new FileEvaluator[AsPairs] {
    override def predictFilesNeededToEvaluate(a: AsPairs, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {

      fileEvaluator.evaluateFilesNeededToEvaluate(a.param, inputs, ioFunctionSet, coerceTo)(fileEvaluator, expressionEvaluator)
    }
  }

  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] = new FileEvaluator[CollectByKey] {
    override def predictFilesNeededToEvaluate(a: CollectByKey, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {

      fileEvaluator.evaluateFilesNeededToEvaluate(a.param, inputs, ioFunctionSet, coerceTo)(fileEvaluator, expressionEvaluator)
    }
  }
}
