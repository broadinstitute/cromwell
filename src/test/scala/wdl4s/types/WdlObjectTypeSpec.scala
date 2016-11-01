package wdl4s.types

import wdl4s.values.{WdlMap, WdlObject, WdlString}
import wdl4s.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Failure, Success}

class WdlObjectTypeSpec extends FlatSpec with Matchers {
  val abcObject = WdlObject(Map(
    "a" -> WdlString("one"),
    "b" -> WdlString("two"),
    "c" -> WdlString("three")
  ))

  val coerceableMap = WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
    WdlString("a") -> WdlString("one"),
    WdlString("b") -> WdlString("two"),
    WdlString("c") -> WdlString("three"))
  )

  val nonCoerceableMap = WdlMap(WdlMapType(WdlStringType, WdlObjectType), Map(
    WdlString("a") -> WdlObject(Map.empty),
    WdlString("b") -> WdlObject(Map.empty),
    WdlString("c") -> WdlObject(Map.empty))
  )
  
  "WdlObject" should "stringify its value" in {
    abcObject.toWdlString shouldEqual "object {a: \"one\", b: \"two\", c: \"three\"}"
  }

  it should "stringify its type" in {
    abcObject.wdlType.toWdlString shouldEqual "Object"
  }

  it should "convert WDL source code to WdlMap" in {
    WdlObjectType.fromWdlString("object {a: \"one\", b: \"two\", c: \"three\"}") shouldEqual abcObject
  }

  it should "coerce a coerceable map into a WdlObject" in {
    WdlObjectType.coerceRawValue(coerceableMap) match {
      case Success(v) => 
        v.wdlType shouldEqual WdlObjectType
        v.toWdlString shouldEqual abcObject.toWdlString
      case Failure(f) => fail("Failed to coerce a map to an object")
    }
  }

  it should "NOT successfully coerce a NON coerceable map into a WdlObject" in {
    WdlObjectType.coerceRawValue(nonCoerceableMap) match {
      case Success(v) => fail("should not have succeeded")
      case Failure(f) => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST" in {
    try {
      WdlObjectType.fromWdlString("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (2)" in {
    try {
      WdlObjectType.fromWdlString("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (3)" in {
    try {
      WdlObjectType.fromWdlString("{1:x(),2:stdout()}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (4)" in {
    try {
      WdlObjectType.fromWdlString("{1:var,2:var}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }


}
