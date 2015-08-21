package cromwell.binding.types

import cromwell.binding.values.{WdlMap, WdlInteger, WdlString}
import cromwell.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsObject, JsArray, JsNumber}

import scala.util.{Failure, Success}

class WdlMapTypeSpec extends FlatSpec with Matchers  {
  val stringIntMap = WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
    WdlString("a") -> WdlInteger(1),
    WdlString("b") -> WdlInteger(2),
    WdlString("c") -> WdlInteger(3)
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
  it should "convert WDL source code to WdlMap" in {
    WdlMapType(WdlStringType, WdlIntegerType).fromWdlString("{\"a\": 1, \"b\": 2, \"c\": 3}") shouldEqual stringIntMap
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWdlString("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (2)" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWdlString("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (3)" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWdlString("{1:x(),2:stdout()}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (4)" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWdlString("{1:var,2:var}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
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
