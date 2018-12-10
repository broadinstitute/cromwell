package wdl.expression

import org.scalatest.{Assertion, FlatSpec, Matchers}
import wdl.draft2.model.expression.WdlStandardLibraryFunctions
import wdl.expression.Draft2SizeFunctionSpec.testFunctions
import wdl.shared.FileSizeLimitationConfig
import wom.expression.EmptyIoFunctionSet
import wom.types._
import wom.values._

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

class Draft2SizeFunctionSpec extends FlatSpec with Matchers {

  behavior of "ReadLikeFunctions.size"

  it should "correctly report a 2048 byte file, in bytes by default" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomSingleFile("blah"))))) { res => assert(res == WomFloat(2048d)) }
  }

  it should "correctly report a 2048 byte file, in bytes" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomSingleFile("blah")), Success(WomString("B"))))) { res => assert(res == WomFloat(2048d)) }
  }

  it should "correctly report a 2048 byte file, in KB" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomSingleFile("blah")), Success(WomString("KB"))))) { res => assert(res == WomFloat(2d)) }
  }

  it should "correctly report a 2048 byte file, in KiB" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomSingleFile("blah")), Success(WomString("Ki"))))) { res => assert(res == WomFloat(2d)) }
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah"))))))) { res => assert(res == WomFloat(2048d)) }
  }

  it should "correctly report the size of a supplied, optional optional, 2048 byte file" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomOptionalValue(
      WomOptionalType(WomSingleFileType),
      Option(WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah"))))
    ))))) { res => assert(res == WomFloat(2048d)) }
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file, in MB" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomOptionalValue(
      WomSingleFileType, Option(WomSingleFile("blah")))), Success(WomString("MB")
    )))) { res => assert(res == WomFloat(0.001953125d)) }
  }

  it should "correctly report that an unsupplied optional file is empty" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomOptionalValue(WomSingleFileType, None))))) { res => assert(res == WomFloat(0d)) }
  }

  it should "correctly report that an unsupplied File?? is empty" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomOptionalValue(WomOptionalType(WomSingleFileType), None))))) { res => assert(res == WomFloat(0d)) }
  }

  it should "correctly report that an unsupplied optional file is empty, even in MB" in {
    val readLike = testFunctions(Success(2048l))
    validate(readLike.size(Seq(Success(WomOptionalValue(WomSingleFileType, None)), Success(WomString("MB")))))  { res => assert(res == WomFloat(0d)) }
  }

  it should "refuse to report file sizes for Ints" in {
    val readLike = testFunctions(Failure(new Exception("Bad result: WdlIntegers shouldn't even be tried for getting file size")))
    val oops = readLike.size(Seq(Success(WomInteger(7))))
    oops match {
      case Success(x) => fail(s"Expected a string to not have a file length but instead got $x")
      case Failure(e) => e.getMessage should be("The 'size' method expects a 'File' or 'File?' argument but instead got Int.")
    }
  }

  it should "refuse to report file sizes for Int?s" in {
    val readLike = testFunctions(Failure(new Exception("Bad result: WdlIntegers shouldn't even be tried for getting file size")))
    val oops = readLike.size(Seq(Success(WomOptionalValue(WomIntegerType, None))))
    oops match {
      case Success(x) => fail(s"Expected a string to not have a file length but instead got $x")
      case Failure(e) => e.getMessage should be("The 'size' method expects a 'File' or 'File?' argument but instead got Int?.")
    }
  }

  it should "pass on underlying size reading errors" in {
    val readLike = testFunctions(Failure(new Exception("'size' inner exception, expect me to be passed on")))
    val oops = readLike.size(Seq(Success(WomSingleFile("blah"))))
    oops match {
      case Success(_) => fail(s"The 'size' engine function didn't return the error generated in the inner 'size' method")
      case Failure(e) => e.getMessage should be("'size' inner exception, expect me to be passed on")
    }
  }

  def validate[A](aTry: Try[A])(f: A => Assertion): Assertion = aTry match {
    case Success(a) => f(a)
    case Failure(e) => fail("Expected success but got failure", e)
  }
}

object Draft2SizeFunctionSpec {
  def testFunctions(sizeResult: Try[Long]): WdlStandardLibraryFunctions = WdlStandardLibraryFunctions.fromIoFunctionSet( new EmptyIoFunctionSet {
    override def size(path: String): Future[Long] = Future.fromTry(sizeResult)
  }, FileSizeLimitationConfig.default)
}
