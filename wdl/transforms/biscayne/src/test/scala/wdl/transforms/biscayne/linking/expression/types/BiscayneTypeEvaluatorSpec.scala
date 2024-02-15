package wdl.transforms.biscayne.linking.expression.types

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.transforms.biscayne.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.biscayne.ast2wdlom._
import wom.types._

class BiscayneTypeEvaluatorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  it should "return nothing from static integer addition" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomIntegerType
    }
  }

  it should "evaluate the right map type from as_map" in {
    val str = """as_map([(1,2), (3,4)])"""
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomMapType(WomIntegerType, WomIntegerType)
    }
  }

  it should "evaluate the right Array[Pair[X, Y]] type from as_pairs" in {
    val str = """as_pairs({ "one": 1, "two": 2, "three": 3 })"""
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomArrayType(WomPairType(WomStringType, WomIntegerType))
    }
  }

  it should "evaluate the type of a sep() function as String" in {
    val str = """ sep(' ', ["a", "b", "c"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomStringType
    }
  }

  it should "evaluate the type of a sep() function with a sub-call to prefix as String" in {
    val str = """ sep(' ', prefix("-i ", ["a", "b", "c"])) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomStringType
    }
  }

  it should "evaluate the type of a suffix() function as Array[String]" in {
    val str = """ suffix('S', ["a", "b", "c"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomArrayType(WomStringType)
    }
  }
}
