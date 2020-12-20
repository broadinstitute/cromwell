package wdl.draft3.transforms.expression

import cats.data.NonEmptyList
import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.linking.expression.types.expressionTypeEvaluator
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types._
import wom.values.WomInteger


class MemberAccessTypeEvaluatorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "member access type evaluator"

  val fiveLiteral = PrimitiveLiteralExpressionElement(WomInteger(5))
  val sixLiteral = PrimitiveLiteralExpressionElement(WomInteger(6))

  it should "find the left and right hand sides of a pair expression" in {
    val pair = PairLiteral(fiveLiteral, sixLiteral)
    val leftExpression: ExpressionElement = ExpressionMemberAccess(pair, NonEmptyList("left", List.empty))
    leftExpression.evaluateType(Map.empty) shouldBeValid WomIntegerType

    val rightExpression: ExpressionElement = ExpressionMemberAccess(pair, NonEmptyList("right", List.empty))
    rightExpression.evaluateType(Map.empty) shouldBeValid WomIntegerType
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
     nestedPairLookup.evaluateType(Map.empty) shouldBeValid WomPairType(WomIntegerType, WomIntegerType)
  }

  it should "evaluate a nested member access on a call output" in {
    val callOutputLookup: ExpressionElement = IdentifierMemberAccess("t", "out", List("left", "right"))
    val linkedValues = Map[UnlinkedConsumedValueHook, GeneratedValueHandle](
      UnlinkedCallOutputOrIdentifierAndMemberAccessHook("t", "out") -> GeneratedCallOutputValueHandle("t", "out", WomPairType(WomPairType(WomIntegerType, WomPairType(WomIntegerType, WomIntegerType)), WomStringType))
    )

    callOutputLookup.evaluateType(linkedValues) shouldBeValid WomPairType(WomIntegerType, WomIntegerType)
  }

  it should "evaluate a nested member access on a struct" in {
    val objectLookup: ExpressionElement = IdentifierMemberAccess("t", "out", List("left", "right"))
    val linkedValues = Map[UnlinkedConsumedValueHook, GeneratedValueHandle](
      UnlinkedCallOutputOrIdentifierAndMemberAccessHook("t", "out") -> GeneratedIdentifierValueHandle(
        linkableName = "t",
        womType = WomCompositeType(Map("out" -> WomPairType(WomPairType(WomIntegerType, WomPairType(WomIntegerType, WomIntegerType)), WomIntegerType)))
      )
    )

    objectLookup.evaluateType(linkedValues) shouldBeValid WomPairType(WomIntegerType, WomIntegerType)
  }

  it should "evaluate a nested member access type on an object" in {
    val objectLookup: ExpressionElement = IdentifierMemberAccess("t", "out", List("left", "right"))
    val linkedValues = Map[UnlinkedConsumedValueHook, GeneratedValueHandle](
      UnlinkedCallOutputOrIdentifierAndMemberAccessHook("t", "out") -> GeneratedIdentifierValueHandle(
        linkableName = "t",
        womType = WomObjectType
      )
    )

    objectLookup.evaluateType(linkedValues) shouldBeValid WomAnyType
  }

  it should "evaluate an identifier lookup" in {
    val identifierLookup: ExpressionElement = IdentifierLookup("foo")
    val linkedValues = Map[UnlinkedConsumedValueHook, GeneratedValueHandle](
      UnlinkedIdentifierHook("foo") -> GeneratedIdentifierValueHandle(linkableName = "foo", womType = WomPairType(WomStringType, WomStringType))
    )

    identifierLookup.evaluateType(linkedValues) shouldBeValid WomPairType(WomStringType, WomStringType)
  }

}
