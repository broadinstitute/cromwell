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
  SQuote,
  SubPosix,
  Suffix,
  Unzip
}
import wdl.model.draft3.graph.expression.FileEvaluator
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.{
  threeParameterFunctionPassthroughFileEvaluator,
  twoParameterFunctionPassthroughFileEvaluator
}
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

object BiscayneFileEvaluators {

  implicit val keysFileEvaluator: FileEvaluator[Keys] = singleParameterPassthroughFileEvaluator
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] = singleParameterPassthroughFileEvaluator
  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] = singleParameterPassthroughFileEvaluator
  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] = singleParameterPassthroughFileEvaluator

  implicit val sepFunctionEvaluator: FileEvaluator[Sep] = twoParameterFunctionPassthroughFileEvaluator[Sep]
  implicit val subPosixFunctionEvaluator: FileEvaluator[SubPosix] =
    threeParameterFunctionPassthroughFileEvaluator[SubPosix]
  implicit val suffixFunctionEvaluator: FileEvaluator[Suffix] = twoParameterFunctionPassthroughFileEvaluator[Suffix]
  implicit val quoteFunctionEvaluator: FileEvaluator[Quote] = singleParameterPassthroughFileEvaluator
  implicit val sQuoteFunctionEvaluator: FileEvaluator[SQuote] = singleParameterPassthroughFileEvaluator

  implicit val minFunctionEvaluator: FileEvaluator[Min] = twoParameterFunctionPassthroughFileEvaluator[Min]
  implicit val maxFunctionEvaluator: FileEvaluator[Max] = twoParameterFunctionPassthroughFileEvaluator[Max]

  implicit val unzipFunctionEvaluator: FileEvaluator[Unzip] =
    singleParameterPassthroughFileEvaluator
}
