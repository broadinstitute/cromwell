package wdl.types

import org.scalatest.{FlatSpec, Matchers}
import wdl.types.WdlFlavoredWomType._
import wdl4s.parser.WdlParser.SyntaxError
import wom.types.{WomArrayType, WomIntegerType}
import wom.values.{WomArray, WomInteger}

class WdlArrayTypeSpec extends FlatSpec with Matchers  {
  val intArray = WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(2), WomInteger(3)))

  it should "convert WDL source code to WdlArray" in {
    WomArrayType(WomIntegerType).fromWorkflowSource("[1,2,3]") shouldEqual intArray
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST" in {
    try {
      WomArrayType(WomIntegerType).fromWorkflowSource("workflow wf{}")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (2)" in {
    try {
      WomArrayType(WomIntegerType).fromWorkflowSource("100")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (3)" in {
    try {
      WomArrayType(WomIntegerType).fromWorkflowSource("[1,stdout()]")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }
  it should "NOT successfully convert WDL source code to WdlArray if passed a bogus AST (4)" in {
    try {
      WomArrayType(WomIntegerType).fromWorkflowSource("[1,var]")
      fail("should not have succeeded")
    } catch {
      case _: SyntaxError => // expected
    }
  }

}
