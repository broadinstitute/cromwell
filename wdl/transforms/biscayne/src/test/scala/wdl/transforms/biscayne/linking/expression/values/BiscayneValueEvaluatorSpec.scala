package wdl.transforms.biscayne.linking.expression.values

import common.assertion.ErrorOrAssertions._
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.transforms.biscayne.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.biscayne.ast2wdlom._
import wom.expression.NoIoFunctionSet
import wom.values.{WomArray, WomInteger, WomMap, WomPair, WomString}
import common.assertion.ManyTimes.intWithTimes

class BiscayneValueEvaluatorSpec extends FlatSpec with Matchers {

  behavior of "biscayne value evaluator"

  // The "did we bring in the base transforms" sanity check:
  it should "return a static integer from static integer addition" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF {
      case e => e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomInteger(6), Seq.empty)
    }
  }

  it should "evaluate an 'as_map' expression correctly" in {
    val str = """ as_map( [("x", 1), ("y", 2), ("z", 3)] ) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedMap: WomMap = WomMap(Map (
      WomString("x") -> WomInteger(1),
      WomString("y") -> WomInteger(2),
      WomString("z") -> WomInteger(3)
    ))

    expr.shouldBeValidPF {
      case e => e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedMap, Seq.empty)
    }
  }

  it should "evaluate an 'as_pairs' expression correctly" in {
    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      val str = """ as_pairs( { 1: "one", 2: "two", 3: three } ) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      val inputs = Map("three" -> WomString("three"))
      val expectedPairs: WomArray = WomArray(Seq(
        WomPair(WomInteger(1), WomString("one")),
        WomPair(WomInteger(2), WomString("two")),
        WomPair(WomInteger(3), WomString("three"))
      ))

      expr.shouldBeValidPF {
        case e => e.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedPairs, Seq.empty)
      }
      ()
    }
  }

  // TODO: Sort out stable map ordering
  // We cannot run this test yet because the map reorders the pairs
  it should "echo correctly via as_pairs(as_map(...))" ignore {

    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      // A value lookup expression:
      val str = """ as_pairs(as_map(echo_me)) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      val expectedPairs: WomArray = WomArray(Seq(
        WomPair(WomInteger(1), WomString("one")),
        WomPair(WomInteger(2), WomString("two")),
        WomPair(WomInteger(3), WomString("three"))
      ))

      val inputs = Map("echo_me" -> expectedPairs)

      expr.shouldBeValidPF {
        case e => e.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedPairs, Seq.empty)
      }
      ()
    }
  }

  it should "fail to evaluate 'as_map' if keys are duplicated" in {

    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      val str = """ as_map( [("x", 1), ("y", 2), ("x", 3)] ) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      expr.shouldBeValidPF {
        case e => e.evaluateValue(Map.empty, NoIoFunctionSet, None).shouldBeInvalid("""Cannot evaluate 'as_map' with duplicated keys: keys can only appear once but "x" appeared 2 times.""")
      }
      ()
    }
  }

  it should "collect by key when keys are duplicated" in {

    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      val str = """ collect_by_key( [("x", 1), ("y", 2), ("x", 3)] ) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      val expectedMap: WomMap = WomMap(Map(
        WomString("x") -> WomArray(Seq(WomInteger(1), WomInteger(3))),
        WomString("y") -> WomArray(Seq(WomInteger(2)))
      ))

      expr.shouldBeValidPF {
        case e => e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedMap, Seq.empty)
      }
      ()
    }
  }

}
