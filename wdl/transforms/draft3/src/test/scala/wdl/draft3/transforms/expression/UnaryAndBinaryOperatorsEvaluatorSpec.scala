package wdl.draft3.transforms.expression

import common.assertion.ErrorOrAssertions._
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wdl.draft3.transforms.linking.expression.types.expressionTypeEvaluator
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.types.{WomBigDecimalType, WomBooleanType, WomIntegerType, WomStringType, WomType}
import wom.values.{WomBigDecimal, WomBoolean, WomFloat, WomInteger, WomString, WomValue}

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

  val womTrue = WomBoolean(true)
  val womFalse = WomBoolean(false)

  val womStringTenCool = WomString("10cool")

  val womBigDecimalTen = WomBigDecimal(10)
  val womBigDecimalEleven = WomBigDecimal(11)
  val womBigDecimalNine = WomBigDecimal(9)
  val womBigDecimalOne = WomBigDecimal(1)
  var womBigDecimalTwentyFive: WomValue = WomBigDecimal(25)
  var womBigDecimalFive: WomValue = WomBigDecimal(5)
  val tenBigDecimalLiteral = PrimitiveLiteralExpressionElement(womBigDecimalTen)
  val nineBigDecimalLiteral = PrimitiveLiteralExpressionElement(womBigDecimalNine)
  val oneBigDecimalLiteral = PrimitiveLiteralExpressionElement(womBigDecimalOne)
  var coolStringLiteral = PrimitiveLiteralExpressionElement(WomString("cool"))
  var twoPointFiveFloatLiteral = PrimitiveLiteralExpressionElement(WomFloat(2.5f))

  // Format:
  // Test name, input expression, output value, output type.
  val expressionTests: List[(String, ExpressionElement, WomValue, WomType)] = List(
    ("+5", UnaryPlus(fiveLiteral), womFive, WomIntegerType),
    ("-5", UnaryNegation(fiveLiteral), WomInteger(-5), WomIntegerType),
    ("-(-5)", UnaryNegation(UnaryNegation(fiveLiteral)), womFive, WomIntegerType),
    ("-(+(-5))", UnaryNegation(UnaryPlus(UnaryNegation(fiveLiteral))), WomInteger(-5), WomIntegerType),
    ("!true", LogicalNot(trueLiteral), womFalse, WomBooleanType),
    ("!(!true)", LogicalNot(LogicalNot(trueLiteral)), womTrue, WomBooleanType),

    ("true || false", LogicalOr(trueLiteral, falseLiteral), womTrue, WomBooleanType),
    ("true && false", LogicalAnd(trueLiteral, falseLiteral), womFalse, WomBooleanType),
    ("true == false", Equals(trueLiteral, falseLiteral), womFalse, WomBooleanType),
    ("true != false", NotEquals(trueLiteral, falseLiteral), womTrue, WomBooleanType),

    ("5 < 6", LessThan(fiveLiteral, sixLiteral), womTrue, WomBooleanType),
    ("5 < 5", LessThan(fiveLiteral, fiveLiteral), womFalse, WomBooleanType),
    ("5 <= 6", LessThanOrEquals(fiveLiteral, sixLiteral), womTrue, WomBooleanType),
    ("5 <= 5", LessThanOrEquals(fiveLiteral, fiveLiteral), womTrue, WomBooleanType),

    ("6 > 5", GreaterThan(sixLiteral, fiveLiteral), womTrue, WomBooleanType),
    ("5 > 5", GreaterThan(fiveLiteral, fiveLiteral), womFalse, WomBooleanType),
    ("6 >= 6", GreaterThanOrEquals(sixLiteral, fiveLiteral), womTrue, WomBooleanType),
    ("5 >= 5", GreaterThanOrEquals(fiveLiteral, fiveLiteral), womTrue, WomBooleanType),

    ("5 + 5", Add(fiveLiteral, fiveLiteral), WomInteger(10), WomIntegerType),
    ("5 - 5", Subtract(fiveLiteral, fiveLiteral), WomInteger(0), WomIntegerType),
    ("5 * 5", Multiply(fiveLiteral, fiveLiteral), WomInteger(25), WomIntegerType),
    ("25 * 5", Divide(PrimitiveLiteralExpressionElement(WomInteger(25)), fiveLiteral), WomInteger(5), WomIntegerType),
    ("27 % 5", Remainder(PrimitiveLiteralExpressionElement(WomInteger(27)), fiveLiteral), WomInteger(2), WomIntegerType),

    //BigDecimal and BigDecimal binary operations
    ("10 + 1", Add(tenBigDecimalLiteral, oneBigDecimalLiteral), womBigDecimalEleven, WomBigDecimalType),
    ("10 * 1", Multiply(tenBigDecimalLiteral, oneBigDecimalLiteral), womBigDecimalTen, WomBigDecimalType),
    ("10 - 1", Subtract(tenBigDecimalLiteral, oneBigDecimalLiteral), womBigDecimalNine, WomBigDecimalType),
    ("10 % 9", Remainder(tenBigDecimalLiteral, nineBigDecimalLiteral), womBigDecimalOne, WomBigDecimalType),
    ("10 / 1", Divide(tenBigDecimalLiteral, oneBigDecimalLiteral), womBigDecimalTen, WomBigDecimalType),
    //Equals, greater and less operations for BigDecimals
    ("10 == 10", Equals(tenBigDecimalLiteral, tenBigDecimalLiteral), womTrue, WomBooleanType),
    ("10 > 1", GreaterThan(tenBigDecimalLiteral, oneBigDecimalLiteral), womTrue, WomBooleanType),
    ("10 < 1", LessThan(tenBigDecimalLiteral, oneBigDecimalLiteral), womFalse, WomBooleanType),
    ("10 >= 1", GreaterThanOrEquals(tenBigDecimalLiteral, oneBigDecimalLiteral), womTrue, WomBooleanType),
    ("10 <= 1", LessThanOrEquals(tenBigDecimalLiteral, oneBigDecimalLiteral), womFalse, WomBooleanType),
    //BigDecimal and Float, Int and String
    ("10 + cool", Add(tenBigDecimalLiteral, coolStringLiteral), womStringTenCool, WomStringType),
    ("10 * 2.5", Multiply(tenBigDecimalLiteral, twoPointFiveFloatLiteral), womBigDecimalTwentyFive, WomBigDecimalType),
    ("10 - 5", Subtract(tenBigDecimalLiteral, fiveLiteral), womBigDecimalFive, WomBigDecimalType)
  )

  expressionTests foreach { case (name, expression, expectedValue, expectedType) =>
    it should s"evaluate the expression '$name'" in {
      expression.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedValue, Seq.empty)
    }

    it should s"evaluate the type of the expression '$name'" in {
      expression.evaluateType(Map.empty) shouldBeValid expectedType
    }
  }
}
