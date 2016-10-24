package wdl4s.types

import wdl4s.values._
import wdl4s.parser.WdlParser.SyntaxError
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import spray.json.JsString

class WdlTypeSpec extends FlatSpec with Matchers {
  "WdlType class" should "stringify WdlBoolean to 'Boolean'" in {
    WdlBooleanType.toWdlString shouldEqual "Boolean"
  }
  it should "stringify WdlInteger to 'Integer'" in {
    WdlIntegerType.toWdlString shouldEqual "Int"
  }
  it should "stringify WdlFloat to 'Float'" in {
    WdlFloatType.toWdlString shouldEqual "Float"
  }
  it should "stringify WdlObject to 'Object'" in {
    WdlObjectType.toWdlString shouldEqual "Object"
  }
  it should "stringify WdlString to 'String'" in {
    WdlStringType.toWdlString shouldEqual "String"
  }
  it should "stringify WdlFile to 'File'" in {
    WdlFileType.toWdlString shouldEqual "File"
  }

  "WdlBoolean" should "support expected coercions" in {
    WdlBooleanType.coerceRawValue("true").get shouldEqual WdlBoolean.True
    WdlBooleanType.coerceRawValue("FALSE").get shouldEqual WdlBoolean.False
    WdlBooleanType.coerceRawValue(false).get shouldEqual WdlBoolean.False

    WdlBooleanType.coerceRawValue("I like turtles").isFailure shouldBe true
  }

  "WdlString" should "support expected coercions" in {
    WdlStringType.coerceRawValue("foo").get shouldEqual WdlString("foo")
    WdlStringType.coerceRawValue(-1).isFailure shouldBe true
  }

  "WdlFile" should "support expected coercions" in {
    WdlFileType.coerceRawValue("/etc/passwd").get shouldEqual WdlFile("/etc/passwd")
    WdlFileType.coerceRawValue(-1).isFailure shouldBe true
  }

  "WdlInteger" should "support expected coercions" in {
    WdlIntegerType.coerceRawValue(42).get shouldEqual WdlInteger(42)
    WdlIntegerType.coerceRawValue("42").get shouldEqual WdlInteger(42)
    WdlIntegerType.coerceRawValue(JsString("42")).get shouldEqual WdlInteger(42)
    WdlIntegerType.coerceRawValue("FAIL").isFailure shouldBe true
  }

  "WdlFloatType" should "support expected coercions" in {
    WdlFloatType.coerceRawValue(33.3).get shouldEqual WdlFloat(33.3)
    WdlFloatType.coerceRawValue("33.3").get shouldEqual WdlFloat(33.3)
    WdlFloatType.coerceRawValue(JsString("33.3")).get shouldEqual WdlFloat(33.3)
    WdlFloatType.coerceRawValue("FAIL").isFailure shouldBe true
  }

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
