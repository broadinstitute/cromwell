package cromwell.binding.types

import cromwell.parser.WdlParser.SyntaxError
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{Matchers, FlatSpec}

class WdlTypeSpec extends FlatSpec with Matchers {
  behavior of "WdlValue"

  val wdlValueRawStrings = Table(
    ("WdlSource", "WdlType"),
    ("String", WdlStringType),
    ("Int", WdlIntegerType),
    ("File", WdlFileType),
    ("Boolean", WdlBooleanType),
    ("Float", WdlFloatType),
    ("Array[Int]", WdlArrayType(WdlIntegerType)),
    ("Array[Array[String]]", WdlArrayType(WdlArrayType(WdlStringType)))
  )

  forAll(wdlValueRawStrings) { (wdlSource, wdlType) =>
    it should s"return the WDL type $wdlType from WDL source string: $wdlSource" in {
      WdlType.fromWdlString(wdlSource) shouldEqual wdlType
    }
  }

  it should "reject an array type without a member type" in {
    try {
      WdlType.fromWdlString("Array")
      fail("Should not be able to parse: `Array`")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "reject an array type with more than one parameterized type" in {
    try {
      WdlType.fromWdlString("Array[String, Int]")
      fail("Should not be able to parse: `Array[String, Int]`")
    } catch {
      case _: SyntaxError => // expected
    }
  }
}
