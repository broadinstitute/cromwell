package wdl.transforms.biscayne.linking.expression.files

import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wdl.transforms.biscayne.Ast2WdlomSpec.{fromString, parser}
import common.assertion.ErrorOrAssertions._
import wom.expression.NoIoFunctionSet
import wom.types.{WomArrayType, WomIntegerType, WomPairType, WomStringType}
import wdl.transforms.biscayne.ast2wdlom._
import wom.values.WomSingleFile
import wdl.transforms.biscayne.linking.expression.values.expressionEvaluator

class BiscayneFileEvaluatorSpec extends FlatSpec with Matchers {
  it should "return nothing from static integer addition" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF {
      case e => e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomIntegerType) shouldBeValid Set.empty
    }
  }

  it should "find the map to read" in {
    val str = """as_pairs(read_map("my_map.txt"))"""
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF {
      case e => e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomArrayType(WomPairType(WomStringType, WomStringType))) shouldBeValid Set(WomSingleFile("my_map.txt"))
    }
  }

  it should "discover the file which would be required to evaluate a sep() function" in {
    val str = """ sep(' ', read_lines("foo.txt")) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF {
      case e => e.predictFilesNeededToEvaluate(Map.empty, NoIoFunctionSet, WomStringType) shouldBeValid Set(WomSingleFile("foo.txt"))
    }
  }
}
