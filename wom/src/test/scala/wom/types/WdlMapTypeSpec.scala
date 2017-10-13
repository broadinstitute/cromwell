package wom.types

import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsNumber, JsObject}
import wom.values.{WdlInteger, WdlMap, WdlObject, WdlString}

import scala.util.{Failure, Success}

class WdlMapTypeSpec extends FlatSpec with Matchers  {
  val stringIntMap = WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
    WdlString("a") -> WdlInteger(1),
    WdlString("b") -> WdlInteger(2),
    WdlString("c") -> WdlInteger(3)
  ))
  val coerceableObject = WdlObject(Map(
    "a" -> WdlString("1"),
    "b" -> WdlString("2"),
    "c" -> WdlString("3")
  ))

  "WdlMap" should "stringify its value" in {
    stringIntMap.toWdlString shouldEqual "{\"a\": 1, \"b\": 2, \"c\": 3}"
  }
  "WdlMapType" should "coerce Map(\"a\":1, \"b\":2, \"c\": 3) into a WdlMap" in {
    WdlMapType(WdlStringType, WdlIntegerType).coerceRawValue(Map("a" -> 1, "b" -> 2, "c" -> 3)) match {
      case Success(array) => array shouldEqual stringIntMap
      case Failure(f) => fail(s"exception while coercing map: $f")
    }
  }
  it should "coerce a JsObject into a WdlMap" in {
    WdlMapType(WdlStringType, WdlIntegerType).coerceRawValue(JsObject(Map("a" -> JsNumber(1), "b" -> JsNumber(2), "c" -> JsNumber(3)))) match {
      case Success(array) => array shouldEqual stringIntMap
      case Failure(f) => fail(s"exception while coercing JsObject: $f")
    }
  }
  it should "stringify its type" in {
    stringIntMap.wdlType.toWdlString shouldEqual "Map[String, Int]"
  }
  it should "coerce a coerceable object into a WdlMap" in {
    WdlMapType(WdlStringType, WdlIntegerType).coerceRawValue(coerceableObject) match {
      case Success(v) =>
        v.wdlType shouldEqual WdlMapType(WdlStringType, WdlIntegerType)
        v.toWdlString shouldEqual stringIntMap.toWdlString
      case Failure(_) => fail("Failed to coerce a map to an object")
    }
  }
  it should "detect invalid map construction if there are mixed types" in {
    try {
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlInteger(0) -> WdlString("foo"), WdlString("x") -> WdlInteger(2)))
      fail("Map initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
  it should "detect invalid map construction if type does not match the input map type" in {
    try {
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlInteger(2) -> WdlInteger(3)))
      fail("Invalid map initialization should have failed")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }
}
