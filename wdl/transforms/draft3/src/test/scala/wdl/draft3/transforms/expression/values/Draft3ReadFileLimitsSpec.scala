package wdl.draft3.transforms.expression.values

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.expression.values.Draft3ReadFileLimitsSpec._
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.transforms.base.linking.expression.values.EngineFunctionEvaluators._
import wom.expression.EmptyIoFunctionSet
import wom.values.WomSingleFile

import scala.concurrent.Future


class Draft3ReadFileLimitsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "ReadLikeFunctions Size Limit Draft 3"
  
  it should "pass correct size limits to the ioFunctions for read_lines" in {
      ReadLines(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
        .evaluateValue(Map.empty, ioFunctionTester(1, ""), None)
        .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_bool" in {
    ReadBoolean(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(2, "true"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_int" in {
    ReadInt(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(3, "0"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_float" in {
    ReadFloat(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(4, "5.0"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_string" in {
    ReadString(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(5, ""), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_json" in {
    ReadJson(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(6, "{}"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_tsv" in {
    ReadTsv(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(7, "a\tb"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_map" in {
    ReadMap(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(8, "a\tb"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_object" in {
    ReadObject(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(9, "a\tb\nc\td"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_objects" in {
    ReadObjects(PrimitiveLiteralExpressionElement(WomSingleFile("blah")))
      .evaluateValue(Map.empty, ioFunctionTester(9, "a\tb\nc\td"), None)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

}

object Draft3ReadFileLimitsSpec {
  def ioFunctionTester(expectedMaxBytes: Int, result: String) = new EmptyIoFunctionSet {
    override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false) = {
      if (maxBytes.contains(expectedMaxBytes)) Future.successful(result)
      else Future.failed(new Exception(s"readFile was called with a max bytes value of ${maxBytes.getOrElse("No value")} but was expecting $expectedMaxBytes"))
    }
  }
}
