package wdl.transforms.biscayne.linking.expression.files

import wdl.model.draft3.elements.ExpressionElement.{
  AsMap,
  AsPairs,
  CollectByKey,
  Keys,
  Max,
  Min,
  Quote,
  Sep,
  SubPosix
  SQuote,
  Suffix
}
import wdl.model.draft3.graph.expression.FileEvaluator
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.{
  threeParameterFunctionPassthroughFileEvaluator,
  twoParameterFunctionPassthroughFileEvaluator
}

object BiscayneFileEvaluators {

  implicit val keysFileEvaluator: FileEvaluator[Keys] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

  implicit val sepFunctionEvaluator: FileEvaluator[Sep] = twoParameterFunctionPassthroughFileEvaluator[Sep]
  implicit val subPosixFunctionEvaluator: FileEvaluator[SubPosix] =
    threeParameterFunctionPassthroughFileEvaluator[SubPosix]
  implicit val suffixFunctionEvaluator: FileEvaluator[Suffix] = twoParameterFunctionPassthroughFileEvaluator[Suffix]
  implicit val quoteFunctionEvaluator: FileEvaluator[Quote] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val sQuoteFunctionEvaluator: FileEvaluator[SQuote] =
    EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

  implicit val minFunctionEvaluator: FileEvaluator[Min] = twoParameterFunctionPassthroughFileEvaluator[Min]
  implicit val maxFunctionEvaluator: FileEvaluator[Max] = twoParameterFunctionPassthroughFileEvaluator[Max]

}
