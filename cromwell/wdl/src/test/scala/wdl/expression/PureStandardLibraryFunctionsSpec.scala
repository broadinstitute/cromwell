package wdl.expression

import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.values._

import scala.util.Success

class PureStandardLibraryFunctionsSpec extends FlatSpec with Matchers {

  behavior of "transpose"

  it should "transpose a 2x3 into a 3x2" in {
    val inArray = WomArray(WomArrayType(WomArrayType(WomIntegerType)), List(
      WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2), WomInteger(3))),
      WomArray(WomArrayType(WomIntegerType), List(WomInteger(4), WomInteger(5), WomInteger(6)))
    ))

    val expectedResult = WomArray(WomArrayType(WomArrayType(WomIntegerType)), List(
      WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(4))),
      WomArray(WomArrayType(WomIntegerType), List(WomInteger(2), WomInteger(5))),
      WomArray(WomArrayType(WomIntegerType), List(WomInteger(3), WomInteger(6)))
    ))

    PureStandardLibraryFunctions.transpose(Seq(Success(inArray))) should be(Success(expectedResult))
  }

  behavior of "length"

  it should "get the right answers" in {

    val two = WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2)))
    PureStandardLibraryFunctions.length(Seq(Success(two))) should be(Success(WomInteger(2)))

    val empty = WomArray(WomArrayType(WomIntegerType), List.empty)
    PureStandardLibraryFunctions.length(Seq(Success(empty))) should be(Success(WomInteger(0)))
  }

  behavior of "prefix"

  it should "prefix things correctly" in {

    val strings = List("foo", "bar", "baz")
    val stringWdlValues = WomArray(WomArrayType(WomStringType), strings map WomString.apply)
    val stringsExpectation = WomArray(WomArrayType(WomStringType), strings map { f => WomString("-f " + f) } )
    PureStandardLibraryFunctions.prefix(Seq(Success(WomString("-f ")), Success(stringWdlValues))) should be(Success(stringsExpectation))

    val noStringWdlValues = WomArray(WomArrayType(WomStringType), List.empty)
    PureStandardLibraryFunctions.prefix(Seq(Success(WomString("-f ")), Success(noStringWdlValues))) should be(Success(WomArray(WomArrayType(WomStringType), Seq.empty)))

    val integers = List(1, 2, 3)
    val integerWdlValues = WomArray(WomArrayType(WomIntegerType), integers map { i => WomInteger.apply(Integer.valueOf(i)) })
    val integersExpectation = WomArray(WomArrayType(WomStringType), integers map { i => WomString("-f " + i)})
    PureStandardLibraryFunctions.prefix(Seq(Success(WomString("-f ")), Success(integerWdlValues))) should be(Success(integersExpectation))
  }

  behavior of "basename"

  List(
    ("my.txt", "my.txt", ".txt", "my"),
    ("/Users/chris/chris.tar.gz", "chris.tar.gz", ".tar.gz", "chris"),
    ("gs://bucket/charlie.bucket", "charlie.bucket", ".wonka", "charlie.bucket")
  ) foreach { case (full, baseWithExtension, suffixToStrip, suffixStripped) =>
    it should s"get the file name for $full" in {
      PureStandardLibraryFunctions.basename(Seq(Success(WomString(full)))) should be(Success(WomString(baseWithExtension)))
    }

    it should s"get the file name for $full and strip the suffix '$suffixToStrip'" in {
      PureStandardLibraryFunctions.basename(Seq(Success(WomString(full)), Success(WomString(suffixToStrip)))) should be(Success(WomString(suffixStripped)))
    }
  }
}
