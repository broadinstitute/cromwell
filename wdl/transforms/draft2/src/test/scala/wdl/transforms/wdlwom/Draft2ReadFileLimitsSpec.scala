package wdl.transforms.wdlwom

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft2.model.{WdlExpression, WdlWomExpression}
import wdl.transforms.wdlwom.Draft2ReadFileLimitsSpec._
import wom.expression.EmptyIoFunctionSet

import scala.concurrent.Future

class Draft2ReadFileLimitsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  behavior of "ReadLikeFunctions Size Limit Draft 2"
  
  it should "pass correct size limits to the ioFunctions for read_lines" in {
    new WdlWomExpression(WdlExpression.fromString("""read_lines("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(1, ""))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_boolean" in {
    new WdlWomExpression(WdlExpression.fromString("""read_boolean("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(2, "true"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_int" in {
    new WdlWomExpression(WdlExpression.fromString("""read_int("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(3, "0"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_float" in {
    new WdlWomExpression(WdlExpression.fromString("""read_float("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(4, "5.0"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_string" in {
    new WdlWomExpression(WdlExpression.fromString("""read_string("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(5, ""))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_json" in {
    new WdlWomExpression(WdlExpression.fromString("""read_json("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(6, "{}"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_tsv" in {
    new WdlWomExpression(WdlExpression.fromString("""read_tsv("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(7, "a\tb"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_map" in {
    new WdlWomExpression(WdlExpression.fromString("""read_map("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(8, "a\tb"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_object" in {
    new WdlWomExpression(WdlExpression.fromString("""read_object("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(9, "a\tb\nc\td"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }

  it should "pass correct size limits to the ioFunctions for read_objects" in {
    new WdlWomExpression(WdlExpression.fromString("""read_objects("blah")"""), null)
      .evaluateValue(Map.empty, ioFunctionTester(9, "a\tb\nc\td"))
      .valueOr(errors => fail(errors.toList.mkString(", ")))
  }
}

object Draft2ReadFileLimitsSpec {
  def ioFunctionTester(expectedMaxBytes: Int, result: String) = new EmptyIoFunctionSet {
    override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false) = {
      if (maxBytes.contains(expectedMaxBytes)) Future.successful(result)
      else Future.failed(new Exception(s"readFile was called with a max bytes value of ${maxBytes.getOrElse("No value")} but was expecting $expectedMaxBytes"))
    }
  }
}
