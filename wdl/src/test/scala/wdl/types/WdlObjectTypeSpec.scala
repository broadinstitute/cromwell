package wdl.types

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.WdlParser.SyntaxError
import wom.types.{WdlMapType, WdlObjectType, WdlStringType}
import wom.values.{WdlMap, WdlObject, WdlString}
import wdl.types.WdlFlavoredWomType._

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

  it should "convert WDL source code to WdlMap" in {
    WdlObjectType.fromWorkflowSource("object {a: \"one\", b: \"two\", c: \"three\"}") shouldEqual abcObject
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST" in {
    try {
      WdlObjectType.fromWorkflowSource("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (2)" in {
    try {
      WdlObjectType.fromWorkflowSource("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (3)" in {
    try {
      WdlObjectType.fromWorkflowSource("{1:x(),2:stdout()}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (4)" in {
    try {
      WdlObjectType.fromWorkflowSource("{1:var,2:var}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

}
