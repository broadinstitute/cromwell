package wdl.transforms.biscayne.linking.expression.files

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{AsMap, AsPairs, CollectByKey, Keys, Max, Min, NoneLiteralElement}
import wdl.model.draft3.graph.expression.{FileEvaluator, ValueEvaluator}
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.twoParameterFunctionPassthroughFileEvaluator
import wom.expression.IoFunctionSet
import wom.types.WomType
import wom.values.{WomFile, WomValue}

object BiscayneFileEvaluators {

  implicit val keysFileEvaluator: FileEvaluator[Keys] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

  implicit val minFunctionEvaluator: FileEvaluator[Min] = twoParameterFunctionPassthroughFileEvaluator[Min]
  implicit val maxFunctionEvaluator: FileEvaluator[Max] = twoParameterFunctionPassthroughFileEvaluator[Max]

  implicit val noneLiteralEvaluator: FileEvaluator[NoneLiteralElement.type] = new FileEvaluator[ExpressionElement.NoneLiteralElement.type] {
    override def predictFilesNeededToEvaluate(a: ExpressionElement.NoneLiteralElement.type,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = Set.empty[WomFile].validNel
  }


}
