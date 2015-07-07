package cromwell.binding.values

import java.nio.file.Paths

import cromwell.binding.{WdlExpression, WdlNamespace}
import cromwell.util.SampleWdl
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

class WdlValueSpec extends FlatSpec with Matchers {

  import TableDrivenPropertyChecks._

  behavior of "WdlValue"

  val wdlValueRawStrings = Table(
    ("wdlValue", "rawString"),
    (WdlBoolean.False, "false"),
    (WdlBoolean.True, "true"),
    (WdlFile("hello/world/path"), "\"hello/world/path\""),
    (WdlFile("hello/world/string"), "\"hello/world/string\""),
    (WdlFloat(0.0), "0.0"),
    (WdlFloat(-0.0), "-0.0"),
    (WdlFloat(Double.PositiveInfinity), "Infinity"),
    (WdlFloat(Double.NegativeInfinity), "-Infinity"),
    (WdlInteger(0), "0"),
    (WdlInteger(Int.MaxValue), "2147483647"),
    (WdlInteger(Int.MinValue), "-2147483648"),
    (WdlString(""), "\"\""),
    (WdlString("test\\'/string"), "\"test\\'/string\""))

  forAll(wdlValueRawStrings) { (wdlValue, rawString) =>
    it should s"exactly convert a ${wdlValue.typeName} to/from WDL source '$rawString'" in {
      val valueAsWdlSource = wdlValue.toWdlString
      valueAsWdlSource should be(rawString)

      val wdlType = wdlValue.wdlType
      val wdlSourceAsValue = wdlType.fromWdlString(valueAsWdlSource)
      wdlSourceAsValue should be(wdlValue)
      wdlSourceAsValue.wdlType should be(wdlType)
    }
  }

  val wdlFloatSpecials = Table(
    ("wdlValue", "rawString", "validateFloat"),
    (WdlFloat(Double.MinPositiveValue), "4.9E-324", { d: Double => d == 0.0 }),
    (WdlFloat(Double.NaN), "NaN", { d: Double => d.isNaN }),
    (WdlFloat(Double.MaxValue), "1.7976931348623157E308", { d: Double => d.isPosInfinity }),
    (WdlFloat(Double.MinValue), "-1.7976931348623157E308", { d: Double => d.isNegInfinity }))

  forAll(wdlFloatSpecials) { (wdlValue, rawString, validateFloat) =>
    it should s"convert a special ${wdlValue.typeName} to/from raw string '$rawString'" in {
      val toRawString = wdlValue.toWdlString
      toRawString should be(rawString)

      val wdlType = wdlValue.wdlType
      val fromRawString = wdlType.fromWdlString(toRawString)
      // Test that this is a special conversion, and is not
      // expected to be equal after a round-trip conversion.
      fromRawString shouldNot be(wdlValue)
      validateFloat(fromRawString.value) should be(right = true)
      fromRawString.wdlType should be(wdlType)
    }
  }

  val wdlExpressionRawStrings = Table(
    ("wdlValue", "rawString"),
    (WdlExpression.fromString(" 1 != 0 "), "1 != 0"),
    (WdlExpression.fromString("10 % 3.5"), "10 % 3.5"),
    (WdlExpression.fromString("10 % 3"), "10 % 3"),
    (WdlExpression.fromString("10-6.7"), "10 - 6.7"),
    (WdlExpression.fromString(""" 1 + "String" """), """1 + "String""""),
    (WdlExpression.fromString("a + b"), "a + b"),
    (WdlExpression.fromString("a(b, c)"), "a(b, c)"),
    (WdlExpression.fromString("\"a\" + \"b\""), "\"a\" + \"b\""),
    (WdlExpression.fromString("a.b.c"), "a.b.c"))

  forAll(wdlExpressionRawStrings) { (wdlValue, rawString) =>
    it should s"resemble a ${wdlValue.typeName} to/from raw string '$rawString'" in {
      val toRawString = wdlValue.toWdlString
      toRawString should be(rawString)

      val wdlType = wdlValue.wdlType
      val fromRawString = wdlType.fromWdlString(toRawString)
      fromRawString shouldNot be(wdlValue)
      fromRawString.toWdlString should be(wdlValue.toWdlString)
      fromRawString.wdlType should be(wdlType)
    }
  }

  val notImplementRawString = Table(
    "wdlValue",
    WdlObject(Map("key" -> WdlString("value"))),
    WdlNamespace.load(SampleWdl.HelloWorld.wdlSource()))

  forAll(notImplementRawString) { wdlValue =>
    it should s"not implement a ${wdlValue.typeName} raw string" in {
      a [NotImplementedError] should be thrownBy wdlValue.toWdlString
      val wdlType = wdlValue.wdlType
      a [NotImplementedError] should be thrownBy wdlType.fromWdlString("")
    }
  }
}
