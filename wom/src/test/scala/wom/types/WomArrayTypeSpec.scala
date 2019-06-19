package wom.types

import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsArray, JsNumber}
import wom.values._

import scala.collection.JavaConverters._
import scala.util.{Failure, Success}

class WomArrayTypeSpec extends FlatSpec with Matchers  {
  val intArray = WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(2), WomInteger(3)))
  "WomArray" should "stringify its value" in {
    intArray.toWomString shouldEqual "[1, 2, 3]"
  }
  it should "tsv serialize an Array[Array[String]]" in {
    val nestedArray = WomArray(
      WomArrayType(WomArrayType(WomStringType)),
      Seq(
        WomArray(WomArrayType(WomStringType), Seq("a", "b").map(WomString)),
        WomArray(WomArrayType(WomStringType), Seq("c", "d").map(WomString))
      )
    )
    nestedArray.tsvSerialize shouldEqual Success("a\tb\nc\td\n")
  }
  it should "fail to TSV serialize an Array[Array[Array[String]]]" in {
    try {
      val array = WomArray(WomArrayType(WomStringType), Seq("a", "b").map(WomString))
      val nestedArray = WomArray(
        WomArrayType(WomArrayType(WomStringType)),
        Seq(
          WomArray(WomArrayType(WomArrayType(WomStringType)), Seq(array, array)),
          WomArray(WomArrayType(WomArrayType(WomStringType)), Seq(array, array))
        )
      )
      nestedArray.tsvSerialize
      fail("should not have succeeded")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  "WomArrayType" should "coerce Seq(1,2,3) into a WomArray" in {
    WomArrayType(WomIntegerType).coerceRawValue(Seq(1,2,3)) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing array: $f")
    }
  }
  it should "coerce a JsArray into a WomArray" in {
    WomArrayType(WomIntegerType).coerceRawValue(JsArray(JsNumber(1), JsNumber(2), JsNumber(3))) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing JsArray: $f")
    }
  }
  it should "coerce a Java List into a WomArray" in {
    WomArrayType(WomIntegerType).coerceRawValue(List(1,2,3).asJava) match {
      case Success(array) => array shouldEqual intArray
      case Failure(f) => fail(s"exception while coercing Java List: $f")
    }
  }
  it should "not coerce single values into one-element arrays" in {
    WomArrayType(WomStringType).coerceRawValue(WomString("edamame is tasty")) match {
      case Success(_) => fail("Unexpected success coercing single value to array")
      case Failure(f) => f.getMessage shouldEqual "No coercion defined from wom value(s) '\"edamame is tasty\"' of type 'String' to 'Array[String]'."
    }
  }
  it should "stringify its type" in {
    intArray.womType.stableName shouldEqual "Array[Int]"
  }
  it should "detect invalid array construction if there are uncoerceable types" in {
    try {
      WomArray(WomArrayType(WomStringType), Seq(WomString("foo"), WomOptionalValue(WomStringType, None)))
      fail("Invalid array initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "detect invalid array construction if type cannot be coerced to the array type" in {
    try {
      WomArray(WomArrayType(WomStringType), Seq(WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(2)))))
      fail("Invalid array initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "correctly coerce values to create an array of optionals from real values" in {
    try {
      val array = WomArray(WomArrayType(WomOptionalType(WomIntegerType)), Seq(WomInteger(2), WomInteger(3)))
      array.value should be(Seq(WomOptionalValue(WomInteger(2)), WomOptionalValue(WomInteger(3))))
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "be able to coerce a map into an array of pairs" in {
    val map = WomMap(WomMapType(WomIntegerType, WomStringType), Map(WomInteger(1) -> WomString("one")))
    val arrayOfPairsType = WomArrayType(WomPairType(WomIntegerType, WomStringType))
    arrayOfPairsType.isCoerceableFrom(map.womType) should be(true)
    arrayOfPairsType.coerceRawValue(map) should be(Success(WomArray(WomArrayType(WomPairType(WomIntegerType, WomStringType)), List(WomPair(WomInteger(1), WomString("one"))))))

  }
}
