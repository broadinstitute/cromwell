package wdl.transforms.biscayne.ast2wdlom

import cats.syntax.validated._
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
    "read_object" -> (_ => "read_object is no longer available in this WDL version. Consider using read_json instead".invalidNel),
    "read_objects" -> (_ => "read_objects is no longer available in this WDL version. Consider using read_json instead".invalidNel),
    "write_object" -> (_ => "write_object is no longer available in this WDL version. Consider using write_json instead".invalidNel),
    "write_objects" -> (_ => "write_objects is no longer available in this WDL version. Consider using write_json instead".invalidNel),
  )
}
