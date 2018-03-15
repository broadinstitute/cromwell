package wdl.draft3.transforms.linking.expression.consumed

import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._

object EngineFunctionEvaluators {
  implicit val rangeFunctionEvaluator: ExpressionValueConsumer[Range] = forOneParamFunction[Range]

  private def forOneParamFunction[A <: OneParamFunctionCallElement] = new ExpressionValueConsumer[A] {
    override def expressionConsumedValueHooks(a: A): Set[UnlinkedConsumedValueHook] = a.param.expressionConsumedValueHooks
  }
}
