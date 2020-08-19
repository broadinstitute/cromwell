package wdl.types


import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft2.model.types.WdlFlavoredWomType._
import wdl.draft2.parser.WdlParser.SyntaxError
import wom.types.{WomIntegerType, WomMapType, WomStringType}
import wom.values.{WomInteger, WomMap, WomString}

class WdlMapTypeSpec extends AnyFlatSpec with Matchers  {
  val stringIntMap = WomMap(WomMapType(WomStringType, WomIntegerType), Map(
    WomString("a") -> WomInteger(1),
    WomString("b") -> WomInteger(2),
    WomString("c") -> WomInteger(3)
  ))
  
  it should "convert WDL source code to WdlMap" in {
    WomMapType(WomStringType, WomIntegerType).fromWorkflowSource("{\"a\": 1, \"b\": 2, \"c\": 3}") shouldEqual stringIntMap
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST" in {
    try {
      WomMapType(WomStringType, WomIntegerType).fromWorkflowSource("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (2)" in {
    try {
      WomMapType(WomStringType, WomIntegerType).fromWorkflowSource("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (3)" in {
    try {
      WomMapType(WomStringType, WomIntegerType).fromWorkflowSource("{1:x(),2:stdout()}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlMap if passed a bogus AST (4)" in {
    try {
      WomMapType(WomStringType, WomIntegerType).fromWorkflowSource("{1:var,2:var}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
}

