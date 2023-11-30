package wdl.transforms.base.linking.expression.files

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}

import scala.util.Try

object UnaryOperatorEvaluators {

  implicit val unaryPlusEvaluator: FileEvaluator[UnaryPlus] = forOperation(_.unaryPlus)
  implicit val unaryNegationEvaluator: FileEvaluator[UnaryNegation] = forOperation(_.unaryMinus)
  implicit val logicalNotEvaluator: FileEvaluator[LogicalNot] = forOperation(_.not)

  private def forOperation[A <: UnaryOperation](op: WomValue => Try[WomValue]) = new FileEvaluator[A] {
    override def predictFilesNeededToEvaluate(a: A,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType
    )(implicit
      fileEvaluator: FileEvaluator[ExpressionElement],
      valueEvaluator: ValueEvaluator[ExpressionElement]
    ): ErrorOr[Set[WomFile]] =
      a.argument.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }
}
