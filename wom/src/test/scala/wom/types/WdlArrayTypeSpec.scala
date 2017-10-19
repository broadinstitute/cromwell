package wom.types

import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsArray, JsNumber}
import wom.values._

import scala.collection.JavaConverters._
import scala.util.{Failure, Success}

class WdlArrayTypeSpec extends FlatSpec with Matchers  {
  val intArray = WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(2), WdlInteger(3)))
  "WdlArray" should "stringify its value" in {
    intArray.toWdlString shouldEqual "[1, 2, 3]"
  }
  it should "tsv serialize an Array[Array[String]]" in {
    val nestedArray = WdlArray(
      WdlArrayType(WdlArrayType(WdlStringType)),
      Seq(
        WdlArray(WdlArrayType(WdlStringType), Seq("a", "b").map(WdlString)),
        WdlArray(WdlArrayType(WdlStringType), Seq("c", "d").map(WdlString))
      )
    )
    nestedArray.tsvSerialize shouldEqual Success("a\tb\nc\td\n")
  }
  it should "fail to TSV serialize an Array[Array[Array[String]]]" in {
    try {
      val array = WdlArray(WdlArrayType(WdlStringType), Seq("a", "b").map(WdlString))
      val nestedArray = WdlArray(
        WdlArrayType(WdlArrayType(WdlStringType)),
        Seq(
          WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), Seq(array, array)),
          WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), Seq(array, array))
        )
      )
      nestedArray.tsvSerialize
      fail("should not have succeeded")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  "WdlArrayType" should "coerce Seq(1,2,3) into a WdlArray" in {
    WdlArrayType(WdlIntegerType).coerceRawValue(Seq(1,2,3)) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing array: $f")
    }
  }
  it should "coerce a JsArray into a WdlArray" in {
    WdlArrayType(WdlIntegerType).coerceRawValue(JsArray(JsNumber(1), JsNumber(2), JsNumber(3))) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing JsArray: $f")
    }
  }
  it should "coerce a Java List into a WdlArray" in {
    WdlArrayType(WdlIntegerType).coerceRawValue(List(1,2,3).asJava) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing Java List: $f")
    }
  }
  it should "coerce single values into one-element arrays" in {
    WdlArrayType(WdlStringType).coerceRawValue(WdlString("edamame is tasty")) match {
      case Success(array) => array shouldEqual WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("edamame is tasty")))
      case Failure(f) => fail("exception coercing single value to array", f)
    }
  }
  it should "stringify its type" in {
    intArray.wdlType.toWdlString shouldEqual "Array[Int]"
  }
  it should "detect invalid array construction if there are uncoerceable types" in {
    try {
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("foo"), WdlOptionalValue(WdlStringType, None)))
      fail("Invalid array initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "detect invalid array construction if type cannot be coerced to the array type" in {
    try {
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(2)))))
      fail("Invalid array initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "correctly coerce values to create an array of optionals from real values" in {
    try {
      val array = WdlArray(WdlArrayType(WdlOptionalType(WdlIntegerType)), Seq(WdlInteger(2), WdlInteger(3)))
      array.value should be(Seq(WdlOptionalValue(WdlInteger(2)), WdlOptionalValue(WdlInteger(3))))
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "be able to coerce a map into an array of pairs" in {
    val map = WdlMap(WdlMapType(WdlIntegerType, WdlStringType), Map(WdlInteger(1) -> WdlString("one")))
    val arrayOfPairsType = WdlArrayType(WdlPairType(WdlIntegerType, WdlStringType))
    arrayOfPairsType.isCoerceableFrom(map.wdlType) should be(true)
    arrayOfPairsType.coerceRawValue(map) should be(Success(WdlArray(WdlArrayType(WdlPairType(WdlIntegerType, WdlStringType)), List(WdlPair(WdlInteger(1), WdlString("one"))))))

  }
}
