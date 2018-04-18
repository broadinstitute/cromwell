package wdl.draft3.transforms.linking.expression.files

import cats.syntax.apply._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.FileEvaluator
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.draft3.transforms.linking.expression.values._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}

object TernaryIfEvaluator {
  implicit val ternaryIfEvaluator: FileEvaluator[TernaryIf] = new FileEvaluator[TernaryIf] {
    override def predictFilesNeededToEvaluate(a: TernaryIf,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType): ErrorOr[Set[WomFile]] = {
      (a.condition.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.ifTrue.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.ifFalse.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ ++ _ }
    }
  }
}
