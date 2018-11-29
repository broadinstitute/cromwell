package wdl.transforms.biscayne.ast2wdlom

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{Keys, AsMap, AsPairs, CollectByKey}
import wdl.transforms.base.ast2wdlom.AstNodeToExpressionElement

object AstToNewExpressionElements {
  val newBiscayneEngineFunctionMakers: Map[String, Vector[ExpressionElement] => ErrorOr[ExpressionElement]] = Map(
    "keys" -> AstNodeToExpressionElement.validateOneParamEngineFunction(Keys, "keys"),
    "as_map" -> AstNodeToExpressionElement.validateOneParamEngineFunction(AsMap, "as_map"),
    "as_pairs" -> AstNodeToExpressionElement.validateOneParamEngineFunction(AsPairs, "as_pairs"),
    "collect_by_key" -> AstNodeToExpressionElement.validateOneParamEngineFunction(CollectByKey, "collect_by_key"),
  )
}
