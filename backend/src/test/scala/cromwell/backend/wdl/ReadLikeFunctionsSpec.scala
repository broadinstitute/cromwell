package cromwell.backend.wdl

import cromwell.core.path.PathBuilder
import org.apache.commons.lang3.NotImplementedException
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.expression.PureStandardLibraryFunctionsLike
import wdl4s.wdl.types.{WdlFileType, WdlIntegerType, WdlOptionalType}
import wdl4s.wdl.values.{WdlFloat, WdlInteger, WdlOptionalValue, WdlSingleFile, WdlString, WdlValue}

import scala.util.{Failure, Success, Try}

class ReadLikeFunctionsSpec extends FlatSpec with Matchers {

  behavior of "ReadLikeFunctions.size"

  it should "correctly report a 2048 byte file, in bytes by default" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlSingleFile("blah")))) should be(Success(WdlFloat(2048d)))
  }

  it should "correctly report a 2048 byte file, in bytes" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlSingleFile("blah")), Success(WdlString("B")))) should be(Success(WdlFloat(2048d)))
  }

  it should "correctly report a 2048 byte file, in KB" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlSingleFile("blah")), Success(WdlString("KB")))) should be(Success(WdlFloat(2.048d)))
  }

  it should "correctly report a 2048 byte file, in KiB" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlSingleFile("blah")), Success(WdlString("Ki")))) should be(Success(WdlFloat(2d)))
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlOptionalValue(WdlFileType, Some(WdlSingleFile("blah")))))) should be(Success(WdlFloat(2048d)))
  }

  it should "correctly report the size of a supplied, optional optional, 2048 byte file" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlOptionalValue(WdlOptionalType(WdlFileType), Some(WdlOptionalValue(WdlFileType, Some(WdlSingleFile("blah")))))))) should be(Success(WdlFloat(2048d)))
  }

  it should "correctly report the size of a supplied, optional, 2048 byte file, in MB" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlOptionalValue(WdlFileType, Some(WdlSingleFile("blah")))), Success(WdlString("MB")))) should be(Success(WdlFloat(0.002048d)))
  }

  it should "correctly report that an unsupplied optional file is empty" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlOptionalValue(WdlFileType, None)))) should be(Success(WdlFloat(0d)))
  }

  it should "correctly report that an unsupplied File?? is empty" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlOptionalValue(WdlOptionalType(WdlFileType), None)))) should be(Success(WdlFloat(0d)))
  }

  it should "correctly report that an unsupplied optional file is empty, even in MB" in {
    val readLike = new TestReadLikeFunctions(Success(2048d))
    readLike.size(Seq(Success(WdlOptionalValue(WdlFileType, None)), Success(WdlString("MB")))) should be(Success(WdlFloat(0d)))
  }

  it should "refuse to report file sizes for Ints" in {
    val readLike = new TestReadLikeFunctions(Failure(new Exception("Bad result: WdlIntegers shouldn't even be tried for getting file size")))
    val oops = readLike.size(Seq(Success(WdlInteger(7))))
    oops match {
      case Success(x) => fail(s"Expected a string to not have a file length but instead got $x")
      case Failure(e) => e.getMessage should be("The 'size' method expects a File argument but instead got Int.")
    }
  }

  it should "refuse to report file sizes for Int?s" in {
    val readLike = new TestReadLikeFunctions(Failure(new Exception("Bad result: WdlIntegers shouldn't even be tried for getting file size")))
    val oops = readLike.size(Seq(Success(WdlOptionalValue(WdlIntegerType, None))))
    oops match {
      case Success(x) => fail(s"Expected a string to not have a file length but instead got $x")
      case Failure(e) => e.getMessage should be("The 'size' method expects a File argument but instead got Int?.")
    }
  }

  it should "pass on underlying size reading errors" in {
    val readLike = new TestReadLikeFunctions(Failure(new Exception("'size' inner exception, expect me to be passed on")))
    val oops = readLike.size(Seq(Success(WdlSingleFile("blah"))))
    oops match {
      case Success(_) => fail(s"The 'size' engine function didn't return the error generated in the inner 'size' method")
      case Failure(e) => e.getMessage should be("'size' inner exception, expect me to be passed on")
    }
  }
}


class TestReadLikeFunctions(sizeResult: Try[Double]) extends PureStandardLibraryFunctionsLike with ReadLikeFunctions {
  override protected def size(file: WdlValue): Try[Double] = sizeResult
  override def pathBuilders: List[PathBuilder] = throw new NotImplementedException("Didn't expect ReadLikefunctionsSpec to need pathBuilders")
}

