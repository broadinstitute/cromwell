package wdl.expression

import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.values._

import scala.util.{Failure, Success}

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

  behavior of "flatten"

  it should "concatenate arrays of arrays of integers" in {
    val ar1 = WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2)))
    val ar2 = WomArray(WomArrayType(WomIntegerType), List(WomInteger(3), WomInteger(4), WomInteger(5)))
    val ar3 = WomArray(WomArrayType(WomIntegerType), List.empty)
    val ar4 = WomArray(WomArrayType(WomIntegerType), List(WomInteger(6)))
    val aar = WomArray(WomArrayType(WomArrayType(WomIntegerType)), List(ar1, ar2, ar3, ar4))
    val flat_ar = WomArray(WomArrayType(WomIntegerType), List(WomInteger(1), WomInteger(2),
                                                              WomInteger(3), WomInteger(4),
                                                              WomInteger(5), WomInteger(6)))
    PureStandardLibraryFunctions.flatten(Seq(Success(aar))) should be(Success(flat_ar))
  }

  it should "concatenate arrays of arrays of strings" in {
    val sar1 = WomArray(WomArrayType(WomStringType), List(WomString("chatting"), WomString("is")))
    val sar2 = WomArray(WomArrayType(WomStringType), List(WomString("great"), WomString("for"), WomString("you")))
    val saar = WomArray(WomArrayType(WomArrayType(WomStringType)), List(sar1, sar2))
    val flat_sar = WomArray(WomArrayType(WomStringType), List(WomString("chatting"), WomString("is"),
                                                              WomString("great"), WomString("for"),
                                                              WomString("you")))
    PureStandardLibraryFunctions.flatten(Seq(Success(saar))) should be(Success(flat_sar))
  }

  it should "return errors for arguments which are arrays with fewer than two dimensions" in {
    val err1 = WomArray(WomArrayType(WomSingleFileType), List.empty)
    PureStandardLibraryFunctions.flatten(Seq(Success(err1))) should be(a[Failure[_]])

    val err2 = WomArray(WomArrayType(WomFloatType), List(WomFloat(1.0), WomFloat(3.2)))
    PureStandardLibraryFunctions.flatten(Seq(Success(err2))) should be(a[Failure[_]])
  }

  it should "return errors for arguments which are not arrays" in {
    val nonArrays: List[WomValue] = List(WomInteger(17),
                                         WomString("banana"),
                                         WomSingleFile("/tmp/bubbles"))
    nonArrays.foreach{ elem =>
      PureStandardLibraryFunctions.flatten(Seq(Success(elem))) should be(a[Failure[_]])
    }
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
