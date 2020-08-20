package wdl.draft3.transforms.expression

import cats.data.NonEmptyList
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.values.{WomInteger, WomObject, WomPair, WomString}


class MemberAccessValueEvaluatorSpec extends AnyFlatSpec with Matchers {

  val fiveLiteral = PrimitiveLiteralExpressionElement(WomInteger(5))
  val sixLiteral = PrimitiveLiteralExpressionElement(WomInteger(6))

  val womFive = WomInteger(5)
  val womSix = WomInteger(6)

  it should "find the left and right hand sides of a pair expression" in {
    val pair = PairLiteral(fiveLiteral, sixLiteral)
    val leftExpression: ExpressionElement = ExpressionMemberAccess(pair, NonEmptyList("left", List.empty))
    leftExpression.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(womFive, Seq.empty)

    val rightExpression: ExpressionElement = ExpressionMemberAccess(pair, NonEmptyList("right", List.empty))
    rightExpression.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(womSix, Seq.empty)
  }

  it should "find the appropriate value in a deeply nested Pair" in {
    val nestedPairLookup: ExpressionElement = ExpressionMemberAccess(
      expression = PairLiteral(
        left = PairLiteral(
          left = fiveLiteral,
          right = PairLiteral(
            left = fiveLiteral,
            right = PairLiteral(
              left = PairLiteral(fiveLiteral, sixLiteral), // This pair is the value we should be selecting
              right = fiveLiteral
            )
          )
        ),
        right = fiveLiteral
      ),
      memberAccessTail = NonEmptyList("left", List("right", "right", "left"))
    )
     nestedPairLookup.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomPair(womFive, womSix), Seq.empty)
  }

  it should "evaluate a nested member access on a call output" in {
    val callOutputLookup: ExpressionElement = IdentifierMemberAccess("t", "out", List("left", "right"))
    val inputs = Map(
      "t.out" -> WomPair(WomPair(womFive, WomPair(womFive, womSix)), womFive)
    )

    callOutputLookup.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomPair(womFive, womSix), Seq.empty)
  }

  it should "evaluate a nested member access on an object" in {
    val objectLookup: ExpressionElement = IdentifierMemberAccess("t", "out", List("left", "right"))
    val inputs = Map(
      "t" -> WomObject(Map("out" -> WomPair(WomPair(womFive, WomPair(womFive, womSix)), womFive)))
    )

    objectLookup.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomPair(womFive, womSix), Seq.empty)
  }

  it should "evaluate an identifier lookup" in {
    val identifierLookup: ExpressionElement = IdentifierLookup("foo")
    val inputs = Map("foo" -> WomString("foo"))

    identifierLookup.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomString("foo"), Seq.empty)
  }

}
