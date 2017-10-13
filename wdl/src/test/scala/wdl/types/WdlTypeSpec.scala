package wdl.types

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.WdlParser.SyntaxError
import wom.types._

class WdlTypeSpec extends FlatSpec with Matchers {

  val wdlValueRawStrings = Table(
    ("WdlSource", "WdlType"),
    ("String", WdlStringType),
    ("Int", WdlIntegerType),
    ("File", WdlFileType),
    ("Boolean", WdlBooleanType),
    ("Float", WdlFloatType),
    ("Array[Int]", WdlArrayType(WdlIntegerType)),
    ("Array[Array[String]]", WdlArrayType(WdlArrayType(WdlStringType))),
    ("Pair[Int, String]", WdlPairType(WdlIntegerType, WdlStringType)),
    ("Pair[Array[Int], String]", WdlPairType(WdlArrayType(WdlIntegerType), WdlStringType))
  )

  forAll(wdlValueRawStrings) { (wdlSource, wdlType) =>
    it should s"return the WDL type $wdlType from WDL source string: $wdlSource" in {
      WdlFlavoredWomType.fromDisplayString(wdlSource) shouldEqual wdlType
    }
  }

  it should "reject an array type without a member type" in {
    try {
      WdlFlavoredWomType.fromDisplayString("Array")
      fail("Should not be able to parse: `Array`")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "reject an array type with more than one parameterized type" in {
    try {
      WdlFlavoredWomType.fromDisplayString("Array[String, Int]")
      fail("Should not be able to parse: `Array[String, Int]`")
    } catch {
      case _: SyntaxError => // expected
    }
  }
}
