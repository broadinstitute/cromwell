package wdl.draft3.transforms.expression

import common.assertion.ErrorOrAssertions._
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.values.{WomBoolean, WomInteger, WomValue}

/**
  * Checks that the draft3 value evaluators for ExpressionElements are wired to forward values through to the appropriate
  * underlying methods on WomValue.
  * ** Not intended as a thorough test of the underlying methods themselves. **
  */
class UnaryAndBinaryOperatorsEvaluatorSpec extends FlatSpec with Matchers{

  val fiveLiteral = PrimitiveLiteralExpressionElement(WomInteger(5))
  val twentyfiveLiteral = PrimitiveLiteralExpressionElement(WomInteger(25))
  val sixLiteral = PrimitiveLiteralExpressionElement(WomInteger(6))
  val trueLiteral = PrimitiveLiteralExpressionElement(WomBoolean(true))
  val falseLiteral = PrimitiveLiteralExpressionElement(WomBoolean(false))

  val womFive = WomInteger(5)
  val womSix = WomInteger(6)

  // Format:
  // Test name, input expression, output value.
  val expressionTests: List[(String, ExpressionElement, WomValue)] = List(
    ("+5", UnaryPlus(fiveLiteral), womFive),
    ("-5", UnaryNegation(fiveLiteral), WomInteger(-5)),
    ("-(-5)", UnaryNegation(UnaryNegation(fiveLiteral)), womFive),
    ("-(+(-5))", UnaryNegation(UnaryPlus(UnaryNegation(fiveLiteral))), WomInteger(-5)),
    ("!true", LogicalNot(trueLiteral), WomBoolean(false)),
    ("!(!true)", LogicalNot(LogicalNot(trueLiteral)), WomBoolean(true)),

    ("true || false", LogicalOr(trueLiteral, falseLiteral), WomBoolean(true)),
    ("true && false", LogicalAnd(trueLiteral, falseLiteral), WomBoolean(false)),
    ("true == false", Equals(trueLiteral, falseLiteral), WomBoolean(false)),
    ("true != false", NotEquals(trueLiteral, falseLiteral), WomBoolean(true)),

    ("5 < 6", LessThan(fiveLiteral, sixLiteral), WomBoolean(true)),
    ("5 < 5", LessThan(fiveLiteral, fiveLiteral), WomBoolean(false)),
    ("5 <= 6", LessThanOrEquals(fiveLiteral, sixLiteral), WomBoolean(true)),
    ("5 <= 5", LessThanOrEquals(fiveLiteral, fiveLiteral), WomBoolean(true)),

    ("6 > 5", GreaterThan(sixLiteral, fiveLiteral), WomBoolean(true)),
    ("5 > 5", GreaterThan(fiveLiteral, fiveLiteral), WomBoolean(false)),
    ("6 >= 6", GreaterThanOrEquals(sixLiteral, fiveLiteral), WomBoolean(true)),
    ("5 >= 5", GreaterThanOrEquals(fiveLiteral, fiveLiteral), WomBoolean(true)),

    ("5 + 5", Add(fiveLiteral, fiveLiteral), WomInteger(10)),
    ("5 - 5", Subtract(fiveLiteral, fiveLiteral), WomInteger(0)),
    ("5 * 5", Multiply(fiveLiteral, fiveLiteral), WomInteger(25)),
    ("25 * 5", Divide(PrimitiveLiteralExpressionElement(WomInteger(25)), fiveLiteral), WomInteger(5)),
    ("27 % 5", Remainder(PrimitiveLiteralExpressionElement(WomInteger(27)), fiveLiteral), WomInteger(2))
  )

  expressionTests foreach { case (name, expression, expected) =>
    it should s"evaluate the expression '$name'" in {
      expression.evaluateValue(Map.empty, NoIoFunctionSet) shouldBeValid expected
    }
  }
}
