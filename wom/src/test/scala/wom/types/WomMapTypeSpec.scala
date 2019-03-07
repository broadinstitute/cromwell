package wom.types

import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsNumber, JsObject}
import wom.values.{WomInteger, WomMap, WomObject, WomString}

import scala.util.{Failure, Success}

class WomMapTypeSpec extends FlatSpec with Matchers  {
  val stringIntMap = WomMap(WomMapType(WomStringType, WomIntegerType), Map(
    WomString("a") -> WomInteger(1),
    WomString("b") -> WomInteger(2),
    WomString("c") -> WomInteger(3)
  ))
  val coerceableObject = WomObject(Map(
    "a" -> WomString("1"),
    "b" -> WomString("2"),
    "c" -> WomString("3")
  ))
  val coerceableTypedObject = WomObject(Map(
    "a" -> WomString("1"),
    "b" -> WomString("2"),
    "c" -> WomString("3")
  ), WomCompositeType(Map("a" -> WomStringType, "b" -> WomStringType, "c" -> WomStringType)))

  "WomMap" should "stringify its value" in {
    stringIntMap.toWomString shouldEqual "{\"a\": 1, \"b\": 2, \"c\": 3}"
  }
  "WomMapType" should "coerce Map(\"a\":1, \"b\":2, \"c\": 3) into a WomMap" in {
    WomMapType(WomStringType, WomIntegerType).coerceRawValue(Map("a" -> 1, "b" -> 2, "c" -> 3)) match {
      case Success(array) => array shouldEqual stringIntMap
      case Failure(f) => fail(s"exception while coercing map: $f")
    }
  }
  it should "coerce a JsObject into a WomMap" in {
    WomMapType(WomStringType, WomIntegerType).coerceRawValue(JsObject(Map("a" -> JsNumber(1), "b" -> JsNumber(2), "c" -> JsNumber(3)))) match {
      case Success(array) => array shouldEqual stringIntMap
      case Failure(f) => fail(s"exception while coercing JsObject: $f")
    }
  }
  it should "stringify its type" in {
    stringIntMap.womType.stableName shouldEqual "Map[String, Int]"
  }
  it should "coerce a coerceable object into a WomMap" in {
    WomMapType(WomStringType, WomIntegerType).coerceRawValue(coerceableObject) match {
      case Success(v) =>
        v.womType shouldEqual WomMapType(WomStringType, WomIntegerType)
        v.toWomString shouldEqual stringIntMap.toWomString
      case Failure(_) => fail("Failed to coerce a map to an object")
    }
  }
  it should "coerce a coerceable typed object into a WomMap" in {
    WomMapType(WomStringType, WomIntegerType).coerceRawValue(coerceableTypedObject) match {
      case Success(v) =>
        v.womType shouldEqual WomMapType(WomStringType, WomIntegerType)
        v.toWomString shouldEqual stringIntMap.toWomString
      case Failure(_) => fail("Failed to coerce a map to an object")
    }
  }

  it should "detect invalid map construction if there are mixed keys" in {
    try {
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomInteger(0) -> WomString("foo"), WomString("x") -> WomString("y")))
      fail("Map initialization should have failed")
    } catch {
      case e: UnsupportedOperationException =>
        e.getMessage shouldEqual "Cannot construct WomMapType(WomStringType,WomStringType) with mixed types in map keys: [WomInteger(0), WomString(x)]"
    }
  }

  it should "detect invalid map construction if there are mixed values" in {
    try {
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("bar") -> WomString("foo"), WomString("x") -> WomInteger(2)))
      fail("Map initialization should have failed")
    } catch {
      case e: UnsupportedOperationException =>
        e.getMessage shouldEqual "Cannot construct WomMapType(WomStringType,WomStringType) with mixed types in map values: [WomString(foo), WomInteger(2)]"
    }
  }

  it should "detect invalid map construction if the keys are all the wrong type" in {
    try {
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomInteger(22) -> WomString("foo"), WomInteger(2222) -> WomString("y")))
      fail("Map initialization should have failed")
    } catch {
      case e: UnsupportedOperationException =>
        e.getMessage shouldEqual "Could not construct a WomMapType(WomStringType,WomStringType) with map keys of unexpected type: [WomInteger(22), WomInteger(2222)]"
    }
  }

  it should "detect invalid map construction if the values are all the wrong type" in {
    try {
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("bar") -> WomInteger(44), WomString("x") -> WomInteger(4444)))
      fail("Map initialization should have failed")
    } catch {
      case e: UnsupportedOperationException =>
        e.getMessage shouldEqual "Could not construct a WomMapType(WomStringType,WomStringType) with map values of unexpected type: [WomInteger(44), WomInteger(4444)]"
    }
  }

  it should "detect invalid map construction if type does not match the input map type" in {
    try {
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomInteger(2) -> WomInteger(3)))
      fail("Invalid map initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
}
