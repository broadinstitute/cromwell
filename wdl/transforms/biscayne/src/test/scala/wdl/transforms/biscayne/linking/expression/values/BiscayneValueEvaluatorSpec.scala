package wdl.transforms.biscayne.linking.expression.values

import common.assertion.ErrorOrAssertions._
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.transforms.biscayne.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.biscayne.ast2wdlom._
import wom.expression.NoIoFunctionSet
import wom.values.{WomArray, WomInteger, WomMap, WomOptionalValue, WomPair, WomString}
import common.assertion.ManyTimes.intWithTimes
import wom.types.{WomIntegerType, WomMapType, WomOptionalType, WomStringType}

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

  it should "evaluate a map literal mixing String?s and Int?s" in {
    val str = """ { "i": i_in, "s": s_in } """
    val inputs = Map(
      "i_in" -> WomOptionalValue(WomIntegerType, Some(WomInteger(1))),
      "s_in" -> WomOptionalValue(WomStringType, Some(WomString("two")))
    )
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedMap: WomMap = WomMap(WomMapType(WomStringType, WomOptionalType(WomStringType)), Map (
      WomString("i") -> WomOptionalValue(WomStringType, Some(WomString("1"))),
      WomString("s") -> WomOptionalValue(WomStringType, Some(WomString("two")))
    ))

    expr.shouldBeValidPF {
      case e => e.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedMap, Seq.empty)
    }
  }

  val escapeTests = Map(
    "\\\\" -> "\\",
    "\\n" -> System.lineSeparator,
    "\\t" -> "\t",
    "\\'" -> "'",
    "\\\"" -> "\"",
    "\\150\\145\\154\\154\\157" -> "hello",
    "\\x68\\x65\\x6C\\x6c\\x6F" -> "hello",
    "\\u0068\\U00000065\\u006C\\U0000006C\\u006F" -> "hello",
    "\\u03A9 (omega)" -> "Ω (omega)"
  )

  escapeTests foreach { case (sequence, expected) =>
    List("\"", "'") foreach { quote =>
      it should s"evaluate the escaping string $quote$sequence$quote as: $expected" in {
        val str = s"$quote$sequence$quote"
        val expectedEvaluation = WomString(expected)
        val expr = fromString[ExpressionElement](str, parser.parse_e)

        expr.shouldBeValidPF {
          case e => e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedEvaluation, Seq.empty)
        }
      }
    }
  }

}
