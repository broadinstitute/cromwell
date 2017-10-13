package wdl.wom_test.types

import org.scalatest.{FlatSpec, Matchers}
import wdl.types.WdlFlavoredWomType._
import wdl4s.parser.WdlParser.SyntaxError
import wom.types.{WdlArrayType, WdlIntegerType}
import wom.values.{WdlArray, WdlInteger}

class WdlArrayTypeSpec extends FlatSpec with Matchers  {
  val intArray = WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(2), WdlInteger(3)))

  it should "convert WDL source code to WdlArray" in {
    WdlArrayType(WdlIntegerType).fromWorkflowSource("[1,2,3]") shouldEqual intArray
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST" in {
    try {
      WdlArrayType(WdlIntegerType).fromWorkflowSource("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (2)" in {
    try {
      WdlArrayType(WdlIntegerType).fromWorkflowSource("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (3)" in {
    try {
      WdlArrayType(WdlIntegerType).fromWorkflowSource("[1,stdout()]")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (4)" in {
    try {
      WdlArrayType(WdlIntegerType).fromWorkflowSource("[1,var]")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

}
