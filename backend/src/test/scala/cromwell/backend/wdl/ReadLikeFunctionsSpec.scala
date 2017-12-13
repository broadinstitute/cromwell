package cromwell.backend.wdl

import cromwell.core.Tags.PostWomTest
import org.scalatest.{FlatSpec, Matchers}
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

class ReadLikeFunctionsSpec extends FlatSpec with Matchers {

  behavior of "ReadLikeFunctions.size"

  it should "correctly report a 2048 byte file, in bytes by default" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomSingleFile("blah")))) should be(Success(WomFloat(2048d)))
  }

  it should "correctly report a 2048 byte file, in bytes" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomSingleFile("blah")), Success(WomString("B")))) should be(Success(WomFloat(2048d)))
  }

  it should "correctly report a 2048 byte file, in KB" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomSingleFile("blah")), Success(WomString("KB")))) should be(Success(WomFloat(2.048d)))
  }

  it should "correctly report a 2048 byte file, in KiB" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomSingleFile("blah")), Success(WomString("Ki")))) should be(Success(WomFloat(2d)))
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah")))))) should
      be(Success(WomFloat(2048d)))
  }

  it should "correctly report the size of a supplied, optional optional, 2048 byte file" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomOptionalValue(
      WomOptionalType(WomSingleFileType),
      Option(WomOptionalValue(WomSingleFileType, Option(WomSingleFile("blah"))))
    )))) should be(Success(WomFloat(2048d)))
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file, in MB" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomOptionalValue(
      WomSingleFileType, Option(WomSingleFile("blah")))), Success(WomString("MB")
    ))) should be(Success(WomFloat(0.002048d)))
  }

  it should "correctly report that an unsupplied optional file is empty" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomOptionalValue(WomSingleFileType, None)))) should be(Success(WomFloat(0d)))
  }

  it should "correctly report that an unsupplied File?? is empty" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomOptionalValue(WomOptionalType(WomSingleFileType), None)))) should
      be(Success(WomFloat(0d)))
  }

  it should "correctly report that an unsupplied optional file is empty, even in MB" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WomOptionalValue(WomSingleFileType, None)), Success(WomString("MB")))) should
      be(Success(WomFloat(0d)))
  }

  it should "refuse to report file sizes for Ints" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Failure(new Exception("Bad result: WdlIntegers shouldn't even be tried for getting file size")))
    val oops = readLike.size(Seq(Success(WomInteger(7))))
    oops match {
      case Success(x) => fail(s"Expected a string to not have a file length but instead got $x")
      case Failure(e) => e.getMessage should be("The 'size' method expects a 'File' or 'File?' argument but instead got Int.")
    }
  }

  it should "refuse to report file sizes for Int?s" taggedAs PostWomTest ignore {
    val readLike = new TestReadLikeFunctions(Failure(new Exception("Bad result: WdlIntegers shouldn't even be tried for getting file size")))
    val oops = readLike.size(Seq(Success(WomOptionalValue(WomIntegerType, None))))
    oops match {
      case Success(x) => fail(s"Expected a string to not have a file length but instead got $x")
      case Failure(e) => e.getMessage should be("The 'size' method expects a 'File' or 'File?' argument but instead got Int?.")
    }
  }

  it should "pass on underlying size reading errors" in {
    val readLike = new TestReadLikeFunctions(Failure(new Exception("'size' inner exception, expect me to be passed on")))
    val oops = readLike.size(Seq(Success(WomSingleFile("blah"))))
    oops match {
      case Success(_) => fail(s"The 'size' engine function didn't return the error generated in the inner 'size' method")
      case Failure(e) => e.getMessage should be("'size' inner exception, expect me to be passed on")
    }
  }
}

// TODO WOM: Fix
class TestReadLikeFunctions(sizeResult: Try[Double]) extends IoFunctionSet {
  override def readFile(path: String): Future[String] = ???

  override def writeFile(path: String, content: String): Future[WomSingleFile] = ???

  override def stdout(params: Seq[Try[WomValue]]): Try[WomSingleFile] = ???

  override def stderr(params: Seq[Try[WomValue]]): Try[WomSingleFile] = ???

  override def glob(pattern: String): Seq[String] = ???

  override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = sizeResult.map(WomFloat.apply)
}
