package wdl.types

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.WdlParser.SyntaxError
import wom.types.{WomMapType, WomObjectType, WomStringType}
import wom.values.{WomMap, WomObject, WomString}
import wdl.types.WdlFlavoredWomType._

class WdlObjectTypeSpec extends FlatSpec with Matchers {
  val abcObject = WomObject(Map(
    "a" -> WomString("one"),
    "b" -> WomString("two"),
    "c" -> WomString("three")
  ))

  val coerceableMap = WomMap(WomMapType(WomStringType, WomStringType), Map(
    WomString("a") -> WomString("one"),
    WomString("b") -> WomString("two"),
    WomString("c") -> WomString("three"))
  )

  val nonCoerceableMap = WomMap(WomMapType(WomStringType, WomObjectType), Map(
    WomString("a") -> WomObject(Map.empty),
    WomString("b") -> WomObject(Map.empty),
    WomString("c") -> WomObject(Map.empty))
  )

  it should "convert WDL source code to WdlMap" in {
    WomObjectType.fromWorkflowSource("object {a: \"one\", b: \"two\", c: \"three\"}") shouldEqual abcObject
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST" in {
    try {
      WomObjectType.fromWorkflowSource("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (2)" in {
    try {
      WomObjectType.fromWorkflowSource("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (3)" in {
    try {
      WomObjectType.fromWorkflowSource("{1:x(),2:stdout()}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (4)" in {
    try {
      WomObjectType.fromWorkflowSource("{1:var,2:var}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

}
