package wdl.transforms.biscayne.linking.expression.files

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.transforms.biscayne.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.biscayne.ast2wdlom._
import wdl.transforms.biscayne.linking.expression.values.expressionEvaluator
import wom.expression.NoIoFunctionSet
import wom.types.{WomArrayType, WomIntegerType, WomPairType, WomStringType}
import wom.values.WomSingleFile

class BiscayneFileEvaluatorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  it should "return nothing from static integer addition" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomIntegerType) shouldBeValid Set.empty
    }
  }

  it should "find the map to read" in {
    val str = """as_pairs(read_map("my_map.txt"))"""
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty,
                                     NoIoFunctionSet,
                                     WomArrayType(WomPairType(WomStringType, WomStringType))
      ) shouldBeValid Set(WomSingleFile("my_map.txt"))
    }
  }

  it should "discover the file which would be required to evaluate a sep() function" in {
    val str = """ sep(' ', read_lines("foo.txt")) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("foo.txt")
      )
    }
  }

  it should "discover the file which would be required to evaluate a sub() function" in {
    val str = """ sub(read_string("my_nice_file.txt"), "foo", "NEW_VAL") """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("my_nice_file.txt")
      )
    }
  }

  it should "discover the file which would be required to evaluate a suffix() function" in {
    val str = """ suffix(' # what a line', read_lines("foo.txt")) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("foo.txt")
      )
    }
  }

  it should "discover the file which would be required to evaluate a quote() function" in {
    val str = """ quote(read_lines("foo.txt")) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("foo.txt")
      )
    }
  }

  it should "discover the file which would be required to evaluate a squote() function" in {
    val str = """ squote(read_lines("foo.txt")) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("foo.txt")
      )
    }
  }

  it should "discover the file which would be required to evaluate a unzip() function" in {
    val str = """ unzip(read_lines("foo.txt")) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("foo.txt")
      )
    }
  }

  it should "discover the file which would be required to evaluate a struct literal" in {
    val str = """ StructWithStringVec{myVec: read_lines("foo.txt")}"""
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(
        WomSingleFile("foo.txt")
      )
    }
  }
}
