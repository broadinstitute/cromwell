package wdl4s.values

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.types.{WdlArrayType, WdlMapType, WdlStringType}
import wdl4s.{SampleWdl, WdlExpression, WdlNamespaceWithWorkflow}

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
    (WdlObject(Map("one" -> WdlString("two"))), "object {one: \"two\"}"),
    (WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("one") -> WdlString("two"))), "{\"one\": \"two\"}"),
    (WdlPair(WdlInteger(1), WdlInteger(2)), "(1, 2)")
  )

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
      validateFloat(fromRawString.value) should be(true)
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

  val testCall = {
    val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.ThreeStep.wdlSource(), Seq.empty).get
    namespace.calls.find(_.unqualifiedName == "wc").get
  }

  val wdlValueMaxedElements = Table(
    ("wdlValue", "maxedElements"),
    (WdlBoolean.False, WdlBoolean.False),
    (WdlFile("hello/world/path"), WdlFile("hello/world/path")),
    (WdlFile("*.txt", isGlob = true), WdlFile("*.txt", isGlob = true)),
    (WdlFloat(0.0), WdlFloat(0.0)),
    (WdlInteger(0), WdlInteger(0)),
    (WdlString(""), WdlString("")),
    (WdlPair(WdlInteger(1), WdlInteger(2)), WdlPair(WdlInteger(1), WdlInteger(2))),
    (WdlOptionalValue(WdlStringType, None), WdlOptionalValue(WdlStringType, None)),
    (
      WdlOptionalValue(WdlStringType, Option(WdlString("optional"))),
      WdlOptionalValue(WdlStringType, Option(WdlString("optional")))
    ),
    (
      WdlObject(Map("0" -> WdlString("zero"))),
      WdlObject(Map("0" -> WdlString("zero")))
    ),
    (
      WdlObject(Map(
        "0" -> WdlString("zero"), "1" -> WdlString("one"), "2" -> WdlString("two"), "3" -> WdlString("three")
      )),
      WdlObject(Map(
        "0" -> WdlString("zero"), "1" -> WdlString("one"), "2" -> WdlString("two")
      ))
    ),
    (
      WdlCallOutputsObject(testCall, Map("0" -> WdlString("zero"))),
      WdlCallOutputsObject(testCall, Map("0" -> WdlString("zero")))
    ),
    (
      WdlCallOutputsObject(testCall, Map(
        "0" -> WdlString("zero"), "1" -> WdlString("one"), "2" -> WdlString("two"), "3" -> WdlString("three")
      )),
      WdlCallOutputsObject(testCall, Map(
        "0" -> WdlString("zero"), "1" -> WdlString("one"), "2" -> WdlString("two")
      ))
    ),
    (
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("0") -> WdlString("zero"))),
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(WdlString("0") -> WdlString("zero")))
    ),
    (
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("0") -> WdlString("zero"),
        WdlString("1") -> WdlString("one"),
        WdlString("2") -> WdlString("two"),
        WdlString("3") -> WdlString("three")
      )),
      WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("0") -> WdlString("zero"),
        WdlString("1") -> WdlString("one"),
        WdlString("2") -> WdlString("two")
      ))
    ),
    (
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("0"))),
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("0")))
    ),
    (
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("0"), WdlString("1"), WdlString("2"), WdlString("3"))),
      WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("0"), WdlString("1"), WdlString("2")))
    ),
    (
      WdlArray(
        WdlArrayType(WdlArrayType(WdlStringType)),
        Seq(
          WdlArray(WdlArrayType(WdlStringType), Seq(
            WdlString("a0"), WdlString("a1"), WdlString("a2"), WdlString("a3")
          )),
          WdlArray(WdlArrayType(WdlStringType), Seq(
            WdlString("b0"), WdlString("b1"), WdlString("b2"), WdlString("b3")
          )),
          WdlArray(WdlArrayType(WdlStringType), Seq(
            WdlString("c0"), WdlString("c1"), WdlString("c2"), WdlString("c3")
          )),
          WdlArray(WdlArrayType(WdlStringType), Seq(
            WdlString("d0"), WdlString("d1"), WdlString("d2"), WdlString("d3")
          ))
        )
      ),
      WdlArray(
        WdlArrayType(WdlArrayType(WdlStringType)),
        Seq(
          WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("a0"), WdlString("a1"), WdlString("a2"))),
          WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("b0"), WdlString("b1"), WdlString("b2"))),
          WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"), WdlString("c1"), WdlString("c2")))
        )
      )
    )
  )

  private def describe(wdlValue: WdlValue): String = {
    wdlValue match {
      case WdlCallOutputsObject(call, outputs) =>
        s"WdlCallOutputsObject(${call.unqualifiedName}, ${outputs.mapValues(_.toWdlString)})"
      case _ => wdlValue.toWdlString
    }
  }

  forAll(wdlValueMaxedElements) { (wdlValue, expected) =>
    it should s"take max elements for ${describe(wdlValue)}" in {
      val actual = WdlValue.takeMaxElements(wdlValue, 3)
      actual should be(expected)
    }
  }
}
