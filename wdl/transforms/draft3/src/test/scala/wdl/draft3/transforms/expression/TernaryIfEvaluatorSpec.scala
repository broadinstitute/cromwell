package wdl.draft3.transforms.expression

import common.assertion.ErrorOrAssertions._
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.transforms.expression.TernaryIfEvaluatorSpec.DontEvaluateMe
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.values.{WomBoolean, WomInteger, WomValue}
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator

/**
  * Checks that the draft3 value evaluators for ExpressionElements are wired to forward values through to the appropriate
  * underlying methods on WomValue.
  * ** Not intended as a thorough test of the underlying methods themselves. **
  */
class TernaryIfEvaluatorSpec extends FlatSpec with Matchers{

  val fiveLiteral = PrimitiveLiteralExpressionElement(WomInteger(5))
  val sixLiteral = PrimitiveLiteralExpressionElement(WomInteger(6))
  val trueLiteral = PrimitiveLiteralExpressionElement(WomBoolean(true))
  val falseLiteral = PrimitiveLiteralExpressionElement(WomBoolean(false))

  val womFive = WomInteger(5)
  val womSix = WomInteger(6)

  // Format:
  // Test name, input expression, output value.
  val expressionTests: List[(String, ExpressionElement, WomValue)] = List(
    ("if true then 5 else ???", TernaryIf(trueLiteral, fiveLiteral, DontEvaluateMe), womFive),
    ("if false then ??? else 6", TernaryIf(falseLiteral, DontEvaluateMe, sixLiteral), womSix)
  )

  expressionTests foreach { case (name, expression, expected) =>
    it should s"evaluate the expression '$name'" in {
      expression.evaluateValue(Map.empty, NoIoFunctionSet) shouldBeValid expected
    }
  }
}

object TernaryIfEvaluatorSpec {

  /**
    * This ExpressionElement won't have any evaluator, so it'll be a failure if we try to evaluate it.
    */
  final case object DontEvaluateMe extends ExpressionElement
}
