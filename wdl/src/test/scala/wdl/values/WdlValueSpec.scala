package wdl.values

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl.{SampleWdl, WdlExpression, WdlNamespaceWithWorkflow}
import wom.types.{WomArrayType, WomMapType, WomStringType}
import wom.values._
import wdl.types.WdlFlavoredWomType._

class WdlValueSpec extends FlatSpec with Matchers {

  import TableDrivenPropertyChecks._

  behavior of "WomValue"

  val wdlValueRawStrings = Table(
    ("womValue", "rawString"),
    (WomBoolean.False, "false"),
    (WomBoolean.True, "true"),
    (WomSingleFile("hello/world/path"), "\"hello/world/path\""),
    (WomSingleFile("hello/world/string"), "\"hello/world/string\""),
    (WomFloat(0.0), "0.0"),
    (WomFloat(-0.0), "-0.0"),
    (WomFloat(Double.PositiveInfinity), "Infinity"),
    (WomFloat(Double.NegativeInfinity), "-Infinity"),
    (WomInteger(0), "0"),
    (WomInteger(Int.MaxValue), "2147483647"),
    (WomInteger(Int.MinValue), "-2147483648"),
    (WomString(""), "\"\""),
    (WomObject(Map("one" -> WomString("two"))), "object {one: \"two\"}"),
    (WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("one") -> WomString("two"))), "{\"one\": \"two\"}"),
    (WomPair(WomInteger(1), WomInteger(2)), "(1, 2)")
  )

  forAll(wdlValueRawStrings) { (womValue, rawString) =>
    it should s"exactly convert a ${womValue.typeName} to/from WDL source '$rawString'" in {
      val valueAsWdlSource = womValue.toWomString
      valueAsWdlSource should be(rawString)

      val womType = womValue.womType
      val wdlSourceAsValue = womType.fromWorkflowSource(valueAsWdlSource)
      wdlSourceAsValue should be(womValue)
      wdlSourceAsValue.womType should be(womType)
    }
  }

  val wdlExpressionRawStrings = Table(
    ("womValue", "rawString"),
    (WdlExpression.fromString(" 1 != 0 "), "1 != 0"),
    (WdlExpression.fromString("10 % 3.5"), "10 % 3.5"),
    (WdlExpression.fromString("10 % 3"), "10 % 3"),
    (WdlExpression.fromString("10-6.7"), "10 - 6.7"),
    (WdlExpression.fromString(""" 1 + "String" """), """1 + "String""""),
    (WdlExpression.fromString("a + b"), "a + b"),
    (WdlExpression.fromString("a(b, c)"), "a(b, c)"),
    (WdlExpression.fromString("\"a\" + \"b\""), "\"a\" + \"b\""),
    (WdlExpression.fromString("a.b.c"), "a.b.c"))

  forAll(wdlExpressionRawStrings) { (womValue, rawString) =>
    it should s"resemble a ${womValue.typeName} to/from raw string '$rawString'" in {
      val toRawString = womValue.toWomString
      toRawString should be(rawString)

      val womType = womValue.womType
      val fromRawString = womType.fromWorkflowSource(toRawString)
      fromRawString shouldNot be(womValue)
      fromRawString.toWomString should be(womValue.toWomString)
      fromRawString.womType should be(womType)
    }
  }

  val testCall = {
    val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.ThreeStep.workflowSource(), Seq.empty).get
    namespace.calls.find(_.unqualifiedName == "wc").get
  }

  val wdlValueMaxedElements = Table(
    ("womValue", "maxedElements"),
    (WomBoolean.False, WomBoolean.False),
    (WomSingleFile("hello/world/path"), WomSingleFile("hello/world/path")),
    (WomGlobFile("*.txt"), WomGlobFile("*.txt")),
    (WomFloat(0.0), WomFloat(0.0)),
    (WomInteger(0), WomInteger(0)),
    (WomString(""), WomString("")),
    (WomPair(WomInteger(1), WomInteger(2)), WomPair(WomInteger(1), WomInteger(2))),
    (WomOptionalValue(WomStringType, None), WomOptionalValue(WomStringType, None)),
    (
      WomOptionalValue(WomStringType, Option(WomString("optional"))),
      WomOptionalValue(WomStringType, Option(WomString("optional")))
    ),
    (
      WomObject(Map("0" -> WomString("zero"))),
      WomObject(Map("0" -> WomString("zero")))
    ),
    (
      WomObject(Map(
        "0" -> WomString("zero"), "1" -> WomString("one"), "2" -> WomString("two"), "3" -> WomString("three")
      )),
      WomObject(Map(
        "0" -> WomString("zero"), "1" -> WomString("one"), "2" -> WomString("two")
      ))
    ),
    (
      WdlCallOutputsObject(testCall, Map("0" -> WomString("zero"))),
      WdlCallOutputsObject(testCall, Map("0" -> WomString("zero")))
    ),
    (
      WdlCallOutputsObject(testCall, Map(
        "0" -> WomString("zero"), "1" -> WomString("one"), "2" -> WomString("two"), "3" -> WomString("three")
      )),
      WdlCallOutputsObject(testCall, Map(
        "0" -> WomString("zero"), "1" -> WomString("one"), "2" -> WomString("two")
      ))
    ),
    (
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("0") -> WomString("zero"))),
      WomMap(WomMapType(WomStringType, WomStringType), Map(WomString("0") -> WomString("zero")))
    ),
    (
      WomMap(WomMapType(WomStringType, WomStringType), Map(
        WomString("0") -> WomString("zero"),
        WomString("1") -> WomString("one"),
        WomString("2") -> WomString("two"),
        WomString("3") -> WomString("three")
      )),
      WomMap(WomMapType(WomStringType, WomStringType), Map(
        WomString("0") -> WomString("zero"),
        WomString("1") -> WomString("one"),
        WomString("2") -> WomString("two")
      ))
    ),
    (
      WomArray(WomArrayType(WomStringType), Seq(WomString("0"))),
      WomArray(WomArrayType(WomStringType), Seq(WomString("0")))
    ),
    (
      WomArray(WomArrayType(WomStringType), Seq(WomString("0"), WomString("1"), WomString("2"), WomString("3"))),
      WomArray(WomArrayType(WomStringType), Seq(WomString("0"), WomString("1"), WomString("2")))
    ),
    (
      WomArray(
        WomArrayType(WomArrayType(WomStringType)),
        Seq(
          WomArray(WomArrayType(WomStringType), Seq(
            WomString("a0"), WomString("a1"), WomString("a2"), WomString("a3")
          )),
          WomArray(WomArrayType(WomStringType), Seq(
            WomString("b0"), WomString("b1"), WomString("b2"), WomString("b3")
          )),
          WomArray(WomArrayType(WomStringType), Seq(
            WomString("c0"), WomString("c1"), WomString("c2"), WomString("c3")
          )),
          WomArray(WomArrayType(WomStringType), Seq(
            WomString("d0"), WomString("d1"), WomString("d2"), WomString("d3")
          ))
        )
      ),
      WomArray(
        WomArrayType(WomArrayType(WomStringType)),
        Seq(
          WomArray(WomArrayType(WomStringType), Seq(WomString("a0"), WomString("a1"), WomString("a2"))),
          WomArray(WomArrayType(WomStringType), Seq(WomString("b0"), WomString("b1"), WomString("b2"))),
          WomArray(WomArrayType(WomStringType), Seq(WomString("c0"), WomString("c1"), WomString("c2")))
        )
      )
    )
  )

  private def describe(womValue: WomValue): String = {
    womValue match {
      case WdlCallOutputsObject(call, outputs) =>
        s"WdlCallOutputsObject(${call.unqualifiedName}, ${outputs.mapValues(_.toWomString)})"
      case _ => womValue.toWomString
    }
  }

  forAll(wdlValueMaxedElements) { (womValue, expected) =>
    it should s"take max elements for ${describe(womValue)}" in {
      val actual = WomValue.takeMaxElements(womValue, 3)
      actual should be(expected)
    }
  }
}
