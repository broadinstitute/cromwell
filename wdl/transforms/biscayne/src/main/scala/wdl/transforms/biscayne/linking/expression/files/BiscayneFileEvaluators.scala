package wdl.transforms.biscayne.linking.expression.files

import wdl.model.draft3.elements.ExpressionElement.{AsMap, AsPairs, CollectByKey, Keys}
import wdl.model.draft3.graph.expression.FileEvaluator
import wdl.transforms.base.linking.expression.files.EngineFunctionEvaluators

object BiscayneFileEvaluators {

  implicit val keysFileEvaluator: FileEvaluator[Keys] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asMapFileEvaluator: FileEvaluator[AsMap] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val asPairsFileEvaluator: FileEvaluator[AsPairs] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator
  implicit val collectByKeyFileEvaluator: FileEvaluator[CollectByKey] = EngineFunctionEvaluators.singleParameterPassthroughFileEvaluator

}
