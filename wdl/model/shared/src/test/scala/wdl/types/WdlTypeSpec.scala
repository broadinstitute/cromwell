package wdl.types

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.WdlParser.SyntaxError
import wom.types._

class WdlTypeSpec extends FlatSpec with Matchers {

  val womValueRawStrings = Table(
    ("DisplayString", "WomType"),
    ("String", WomStringType),
    ("Int", WomIntegerType),
    ("File", WomSingleFileType),
    ("Boolean", WomBooleanType),
    ("Float", WomFloatType),
    ("Array[Int]", WomArrayType(WomIntegerType)),
    ("Array[Array[String]]", WomArrayType(WomArrayType(WomStringType))),
    ("Pair[Int, String]", WomPairType(WomIntegerType, WomStringType)),
    ("Pair[Array[Int], String]", WomPairType(WomArrayType(WomIntegerType), WomStringType))
  )

  forAll(womValueRawStrings) { (displayString, womType) =>
    it should s"return the WDL type $womType from WOM type display string: $displayString" in {
      WdlFlavoredWomType.fromDisplayString(displayString) shouldEqual womType
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
