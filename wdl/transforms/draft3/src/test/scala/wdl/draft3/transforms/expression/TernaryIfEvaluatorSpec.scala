package wdl.draft3.transforms.expression

import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.expression.TernaryIfEvaluatorSpec.DontEvaluateMe
import wdl.draft3.transforms.linking.expression.types.expressionTypeEvaluator
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.types.{WomFloatType, WomIntegerType, WomType}
import wom.values.{WomBoolean, WomFloat, WomInteger, WomValue}


/**
  * Checks that the draft3 value evaluators for ExpressionElements are wired to forward values through to the appropriate
  * underlying methods on WomValue.
  * ** Not intended as a thorough test of the underlying methods themselves. **
  */
class TernaryIfEvaluatorSpec extends AnyFlatSpec with Matchers {

  val fiveLiteral = PrimitiveLiteralExpressionElement(WomInteger(5))
  val sixLiteral = PrimitiveLiteralExpressionElement(WomInteger(6))
  val trueLiteral = PrimitiveLiteralExpressionElement(WomBoolean(true))
  val falseLiteral = PrimitiveLiteralExpressionElement(WomBoolean(false))

  val womFive = WomInteger(5)
  val womSix = WomInteger(6)

  // Format:
  // Test name, input expression, output value.
  val valueTests: List[(String, ExpressionElement, WomValue)] = List(
    ("if true then 5 else ???", TernaryIf(trueLiteral, fiveLiteral, DontEvaluateMe), womFive),
    ("if false then ??? else 6", TernaryIf(falseLiteral, DontEvaluateMe, sixLiteral), womSix)
  )

  // Format:
  // Test name, input expression, output type.
  val typeTests: List[(String, ExpressionElement, WomType)] = List(
    ("if true then 5 else 6", TernaryIf(trueLiteral, fiveLiteral, sixLiteral), WomIntegerType),
    ("if false then 5.5 else 6", TernaryIf(falseLiteral, PrimitiveLiteralExpressionElement(WomFloat(5.5)), sixLiteral), WomFloatType)
  )

  valueTests foreach { case (name, expression, expected) =>
    it should s"evaluate the expression '$name'" in {
      expression.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expected, Seq.empty)
    }
  }

  typeTests foreach { case (name, expression, expected) =>
    it should s"evaluate the expression '$name'" in {
      expression.evaluateType(Map.empty) shouldBeValid expected
    }
  }
}

object TernaryIfEvaluatorSpec {

  /**
    * This ExpressionElement won't have any evaluator, so it'll be a failure if we try to evaluate it.
    */
  final case object DontEvaluateMe extends ExpressionElement
}
