package wdl4s.expression

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlInteger, WdlString}

import scala.util.Success

class PureStandardLibraryFunctionsSpec extends FlatSpec with Matchers {

  behavior of "transpose"

  it should "transpose a 2x3 into a 3x2" in {
    val inArray = WdlArray(WdlArrayType(WdlArrayType(WdlIntegerType)), List(
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(2), WdlInteger(3))),
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(4), WdlInteger(5), WdlInteger(6)))
    ))

    val expectedResult = WdlArray(WdlArrayType(WdlArrayType(WdlIntegerType)), List(
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(4))),
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(2), WdlInteger(5))),
      WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(3), WdlInteger(6)))
    ))

    PureStandardLibraryFunctions.transpose(Seq(Success(inArray))) should be(Success(expectedResult))
  }

  behavior of "length"

  it should "get the right answers" in {

    val two = WdlArray(WdlArrayType(WdlIntegerType), List(WdlInteger(1), WdlInteger(2)))
    PureStandardLibraryFunctions.length(Seq(Success(two))) should be(Success(WdlInteger(2)))

    val empty = WdlArray(WdlArrayType(WdlIntegerType), List.empty)
    PureStandardLibraryFunctions.length(Seq(Success(empty))) should be(Success(WdlInteger(0)))
  }

  behavior of "prefix"

  it should "prefix things correctly" in {

    val strings = List("foo", "bar", "baz")
    val stringWdlValues = WdlArray(WdlArrayType(WdlStringType), strings map WdlString.apply)
    val stringsExpectation = WdlArray(WdlArrayType(WdlStringType), strings map { f => WdlString("-f " + f) } )
    PureStandardLibraryFunctions.prefix(Seq(Success(WdlString("-f ")), Success(stringWdlValues))) should be(Success(stringsExpectation))

    val noStringWdlValues = WdlArray(WdlArrayType(WdlStringType), List.empty)
    PureStandardLibraryFunctions.prefix(Seq(Success(WdlString("-f ")), Success(noStringWdlValues))) should be(Success(WdlArray(WdlArrayType(WdlStringType), Seq.empty)))

    val integers = List(1, 2, 3)
    val integerWdlValues = WdlArray(WdlArrayType(WdlIntegerType), integers map { i => WdlInteger.apply(Integer.valueOf(i)) })
    val integersExpectation = WdlArray(WdlArrayType(WdlStringType), integers map { i => WdlString("-f " + i)})
    PureStandardLibraryFunctions.prefix(Seq(Success(WdlString("-f ")), Success(integerWdlValues))) should be(Success(integersExpectation))
  }

  behavior of "basename"

  List(
    ("my.txt", "my.txt", ".txt", "my"),
    ("/Users/chris/chris.tar.gz", "chris.tar.gz", ".tar.gz", "chris"),
    ("gs://bucket/charlie.bucket", "charlie.bucket", ".wonka", "charlie.bucket")
  ) foreach { case (full, baseWithExtension, suffixToStrip, suffixStripped) =>
    it should s"get the file name for $full" in {
      PureStandardLibraryFunctions.basename(Seq(Success(WdlString(full)))) should be(Success(WdlString(baseWithExtension)))
    }

    it should s"get the file name for $full and strip the suffix '$suffixToStrip'" in {
      PureStandardLibraryFunctions.basename(Seq(Success(WdlString(full)), Success(WdlString(suffixToStrip)))) should be(Success(WdlString(suffixStripped)))
    }
  }
}
