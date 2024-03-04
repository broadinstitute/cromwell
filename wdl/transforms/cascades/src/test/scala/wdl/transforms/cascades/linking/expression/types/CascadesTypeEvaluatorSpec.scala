package wdl.transforms.cascades.linking.expression.types

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.transforms.cascades.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.cascades.ast2wdlom._
import wom.types._

class CascadesTypeEvaluatorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
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

  it should "evaluate the type of a sub() function as String" in {
    val str = """ sub("input", "^pattern$", "s") """
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

  it should "evaluate the type of a quote() function as Array[String]" in {
    val str = """ quote([1, 2, 3]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomArrayType(WomStringType)
    }
  }

  it should "evaluate the type of a quote() function with an empty array as Array[Nothing]" in {
    val str = """ quote([]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomArrayType(WomNothingType)
    }
  }

  it should "evaluate the type of a squote() function as Array[String]" in {
    val str = """ squote([1, 2, 3]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomArrayType(WomStringType)
    }
  }

  it should "evaluate the type of an squote() function with an empty array as Array[Nothing]" in {
    val str = """ squote([]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomArrayType(WomNothingType)
    }
  }

  it should "evaluate the type of an unzip() function as Pair[Array[X], Array[Y]]" in {
    val string_and_int = """ unzip([("one", 1),("two", 2),("three", 3)]) """
    val string_and_int_expr = fromString[ExpressionElement](string_and_int, parser.parse_e)
    string_and_int_expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomPairType(WomArrayType(WomStringType), WomArrayType(WomIntegerType))
    }

    val int_and_int = """ unzip([(1,2),(3,4),(5,6)]) """
    val int_and_int_expr = fromString[ExpressionElement](int_and_int, parser.parse_e)
    int_and_int_expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomPairType(WomArrayType(WomIntegerType), WomArrayType(WomIntegerType))
    }
  }

  it should "evaluate the type of an unzip() function on an empty collection as Pair[Array[Any], Array[Any]]" in {
    val empty = """ unzip([]) """
    val empty_unzip_expr = fromString[ExpressionElement](empty, parser.parse_e)
    empty_unzip_expr.shouldBeValidPF { case e =>
      e.evaluateType(Map.empty) shouldBeValid WomPairType(WomArrayType(WomAnyType), WomArrayType(WomAnyType))
    }
  }
}
