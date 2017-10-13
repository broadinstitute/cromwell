package wdl.types


import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.WdlParser.SyntaxError
import wom.types.{WdlIntegerType, WdlMapType, WdlStringType}
import wom.values.{WdlInteger, WdlMap, WdlString}
import wdl.types.WdlFlavoredWomType._

class WdlMapTypeSpec extends FlatSpec with Matchers  {
  val stringIntMap = WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
    WdlString("a") -> WdlInteger(1),
    WdlString("b") -> WdlInteger(2),
    WdlString("c") -> WdlInteger(3)
  ))
  
  it should "convert WDL source code to WdlMap" in {
    WdlMapType(WdlStringType, WdlIntegerType).fromWorkflowSource("{\"a\": 1, \"b\": 2, \"c\": 3}") shouldEqual stringIntMap
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWorkflowSource("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (2)" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWorkflowSource("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (3)" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWorkflowSource("{1:x(),2:stdout()}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (4)" in {
    try {
      WdlMapType(WdlStringType, WdlIntegerType).fromWorkflowSource("{1:var,2:var}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
}

