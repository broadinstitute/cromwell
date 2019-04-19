package wdl.draft3.transforms.expression.values

import cats.data.Validated.{Invalid, Valid}
import common.validation.ErrorOr.ErrorOr
import org.scalatest.{Assertion, FlatSpec, Matchers}
import wdl.draft3.transforms.expression.values.Draft3SizeFunctionSpec.testFunctions
import wdl.transforms.base.linking.expression.values.EngineFunctionEvaluators.sizeFunctionEvaluator
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement, Size}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.draft3.transforms.linking.expression.values.expressionEvaluator
import wom.expression.{EmptyIoFunctionSet, IoFunctionSet}
import wom.types._
import wom.values._

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

class Draft3SizeFunctionSpec extends FlatSpec with Matchers {

  behavior of "ReadLikeFunctions.size"

  it should "correctly report a 2048 byte file, in bytes by default" in {
    validate(Size(PrimitiveLiteralExpressionElement(WomSingleFile("blah")), None).evaluateValue(Map.empty, testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2048d))
    }
  }

  it should "correctly report a 2048 byte file, in bytes" in {
    validate(Size(PrimitiveLiteralExpressionElement(WomSingleFile("blah")), Some(PrimitiveLiteralExpressionElement(WomString("B")))).evaluateValue(Map.empty, testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2048d))
    }
  }

  it should "correctly report a 2048 byte file, in KB" in {
    validate(Size(PrimitiveLiteralExpressionElement(WomSingleFile("blah")), Some(PrimitiveLiteralExpressionElement(WomString("KB")))).evaluateValue(Map.empty, testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2.0d))
    }
  }

  it should "correctly report a 2048 byte file, in KiB" in {
    validate(Size(PrimitiveLiteralExpressionElement(WomSingleFile("blah")), Some(PrimitiveLiteralExpressionElement(WomString("KiB")))).evaluateValue(Map.empty, testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2d))
    }
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file" in {
    val value = WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah")))

    validate(Size(IdentifierLookup("x"), None).evaluateValue(Map("x" -> value), testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2048d))
    }
  }

  it should "correctly report the size of a supplied, optional optional, 2048 byte file" in {
    val value = WomOptionalValue(WomOptionalType(WomSingleFileType), Option(WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah")))))

    validate(Size(IdentifierLookup("x"), None).evaluateValue(Map("x" -> value), testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2048d))
    }
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file, in MB" in {
    val value = WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah")))

    validate(Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("MB")))).evaluateValue(Map("x" -> value), testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(0.001953125d))
    }
  }

  it should "correctly report that an unsupplied optional file is empty" in {
    val value = WomOptionalValue(WomSingleFileType, None)

    validate(Size(IdentifierLookup("x"), None).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("Bad call to size on an empty optional"))), None)) {
      res => assert(res.value == WomFloat(0d))
    }
  }

  it should "correctly report that an unsupplied File?? is empty" in {
      val value = WomOptionalValue(WomOptionalType(WomSingleFileType), Option(WomOptionalValue(WomSingleFileType, None)))

      validate(Size(IdentifierLookup("x"), None).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("Bad call to size on an empty optional"))), None)) {
        res => assert(res.value == WomFloat(0d))
      }
  }

  it should "correctly report that an unsupplied optional file is empty, even in MB" in {
    val value = WomOptionalValue(WomSingleFileType, None)

    validate(Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("MB")))).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("Bad call to size on an empty optional"))), None)) {
      res => assert(res.value == WomFloat(0d))
    }
  }

  it should "correctly report the size of an array of files, in GiB" in {
    val value = WomArray(Seq(WomSingleFile("blah"), WomSingleFile("blah")))

    validate(Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("GiB")))).evaluateValue(Map("x" -> value), testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2048d * 2 / 1024 / 1024 / 1024))
    }
  }

  it should "correctly report that the size of an array of unsupplied optional files is empty, in MB" in {
    val value = WomArray(Seq(WomOptionalValue(WomSingleFileType, None), WomOptionalValue(WomSingleFileType, None)))

    validate(Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("MB")))).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("Bad call to size on an empty optional"))), None)) {
      res => assert(res.value == WomFloat(0d))
    }
  }

  it should "correctly report the size of a mixed Array[File?] - some supplied and some not" in {
    val value = WomArray(Seq(WomOptionalValue(WomSingleFileType, Some(WomSingleFile("blah"))), WomOptionalValue(WomSingleFileType, None), WomOptionalValue(WomSingleFileType, Some(WomSingleFile("blah"))), WomOptionalValue(WomSingleFileType, None)))

    validate(Size(IdentifierLookup("x"), None).evaluateValue(Map("x" -> value), testFunctions(Success(2048l)), None)) {
      res => assert(res.value == WomFloat(2048d * 2))
    }
  }

  it should "refuse to report file sizes for Ints" in {
    val value = WomInteger(55)
    val oops = Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("MB")))).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("Bad call to size on an Int"))), None)
    oops match {
      case Valid(x) => fail(s"Expected an Integer to not have a file length but instead got ${x.value.toWomString}")
      case Invalid(e) => e.head should be("The 'size' method expects a 'File', 'File?', 'Array[File]' or Array[File?] argument but instead got Int.")
    }
  }

  it should "refuse to report file sizes for Int?s" in {
    val value = WomOptionalValue(WomIntegerType, Some(WomInteger(55)))
    val oops = Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("MB")))).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("Bad call to size on an Int?"))), None)
    oops match {
      case Valid(x) => fail(s"Expected an Int? to not have a file length but instead got ${x.value.toWomString}")
      case Invalid(e) => e.head should be("The 'size' method expects a 'File', 'File?', 'Array[File]' or Array[File?] argument but instead got Int?.")
    }
  }

  it should "pass on underlying size reading errors" in {
    val value = WomOptionalValue(WomSingleFileType, Some(WomSingleFile("blah")))
    val oops = Size(IdentifierLookup("x"), Some(PrimitiveLiteralExpressionElement(WomString("MB")))).evaluateValue(Map("x" -> value), testFunctions(Failure(new Exception("'size' inner exception, expect me to be passed on"))), None)
    oops match {
      case Valid(x) => fail(s"Expected an Int? to not have a file length but instead got ${x.value.toWomString}")
      case Invalid(e) => e.head should be("'size' inner exception, expect me to be passed on")
    }
  }

  def validate[A](errorOr: ErrorOr[A])(f: A => Assertion): Assertion = errorOr match {
    case Valid(value) => f(value)
    case Invalid(e) => fail(s"Expected success but got failures: [${e.toList.mkString(", ")}]")
  }
}

object Draft3SizeFunctionSpec {
  def testFunctions(sizeResult: Try[Long]): IoFunctionSet = new EmptyIoFunctionSet {
    override def size(path: String): Future[Long] = Future.fromTry(sizeResult)
  }
}
