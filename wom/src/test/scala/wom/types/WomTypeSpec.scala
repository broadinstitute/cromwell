package wom.types

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import spray.json.JsString
import wom.values._

import scala.runtime.ScalaRunTime
import scala.util.Random

class WomTypeSpec extends FlatSpec with Matchers {
  "WomType class" should "stringify WomBoolean to 'Boolean'" in {
    WomBooleanType.stableName shouldEqual "Boolean"
  }
  it should "stringify WomInteger to 'Integer'" in {
    WomIntegerType.stableName shouldEqual "Int"
  }
  it should "stringify WomFloat to 'Float'" in {
    WomFloatType.stableName shouldEqual "Float"
  }
  it should "stringify WomObject to 'Object'" in {
    WomObjectType.stableName shouldEqual "Object"
  }
  it should "stringify WomString to 'String'" in {
    WomStringType.stableName shouldEqual "String"
  }
  it should "stringify WomFile to 'File'" in {
    WomSingleFileType.stableName shouldEqual "File"
  }

  val rawValuesCoercedToType = Table(
    (
      "Raw Value",
      "WomType",
      "Exception Class",
      "Exception Message Regex"
    ),
    (
      WomString("hello"),
      WomIntegerType,
      classOf[NumberFormatException],
      "For input string: \"hello\""
    ),
    (
      WomInteger(0),
      WomBooleanType,
      classOf[IllegalArgumentException],
      """No coercion defined from wom value\(s\) '0' of type 'Int' to 'Boolean'."""
    ),
    (
      0,
      WomBooleanType,
      classOf[IllegalArgumentException],
      "No coercion defined from '0' of type 'java.lang.Integer' to 'Boolean'."
    ),
    (
      Array(0, 1, 2, 3, 4),
      WomBooleanType,
      classOf[IllegalArgumentException],
      """No coercion defined from 'Array\(0, 1, 2\)' of type 'int\[\]' to 'Boolean'."""
    ),
    (
      new AnyRef {},
      WomBooleanType,
      classOf[IllegalArgumentException],
      "No coercion defined from" +
        """ 'wom.types.WomTypeSpec\$\$anon\$(.*)@.*' of type""" +
        """ 'wom.types.WomTypeSpec\$\$anon\$\1' to 'Boolean'."""
    ),
    (
      WomArray(WomArrayType(WomOptionalType(WomIntegerType)), Seq(
        WomOptionalValue(WomInteger(0)),
        WomOptionalValue(WomInteger(1)),
        WomOptionalValue(WomInteger(2)),
        WomOptionalValue(WomInteger(3)),
        WomOptionalValue(WomInteger(4)))
      ),
      WomOptionalType(WomMaybeEmptyArrayType(WomIntegerType)),
      classOf[IllegalArgumentException],
      """No coercion defined from wom value\(s\) '\[0, 1, 2\]' of type 'Array\[Int\?\]' to 'Array\[Int\]\?'."""
    ),
    (
      WomArray(WomArrayType(WomOptionalType(WomIntegerType)), Seq(WomOptionalValue.none(WomIntegerType))),
      WomOptionalType(WomMaybeEmptyArrayType(WomIntegerType)),
      classOf[IllegalArgumentException],
      """No coercion defined from wom value\(s\) '\[null\]' of type 'Array\[Int\?\]' to 'Array\[Int\]\?'."""
    ),
    (
      WomOptionalValue.none(WomArrayType(WomIntegerType)),
      WomMaybeEmptyArrayType(WomOptionalType(WomIntegerType)),
      classOf[IllegalArgumentException],
      """No coercion defined from wom value\(s\) 'null' of type 'Array\[Int\]\?' to 'Array\[Int\?\]'."""
    )
  )

  private def describe(any: Any): String = {
    any match {
      case womValue: WomValue => s"wom value ${womValue.toWomString} of type ${womValue.womType.stableName}"
      case _ => s"scala value ${ScalaRunTime.stringOf(any)}"
    }
  }

  forAll(rawValuesCoercedToType) { (rawValue, womType, exceptionClass, exceptionMessage) =>
    it should s"fail coercing ${womType.stableName} from ${describe(rawValue)}" in {
      val exception = womType.coerceRawValue(rawValue).failed.get
      exception.getClass should be(exceptionClass)
      exception.getMessage should fullyMatch regex exceptionMessage
    }
  }

  "WomBoolean" should "support expected coercions" in {
    WomBooleanType.coerceRawValue("true").get shouldEqual WomBoolean.True
    WomBooleanType.coerceRawValue("FALSE").get shouldEqual WomBoolean.False
    WomBooleanType.coerceRawValue(false).get shouldEqual WomBoolean.False
    WomBooleanType.coerceRawValue(WomOptionalValue(WomBooleanType, Option(WomBoolean(true)))).get shouldEqual WomBoolean.True
    WomBooleanType.coerceRawValue("I like turtles").isFailure shouldBe true
  }

  "WomString" should "support expected coercions" in {
    WomStringType.coerceRawValue("foo").get shouldEqual WomString("foo")
    WomStringType.coerceRawValue(WomOptionalValue(WomStringType, Option(WomString("foo")))).get shouldEqual WomString("foo")
    WomStringType.coerceRawValue(-1).isFailure shouldBe true
  }

  "WomFile" should "support expected coercions" in {
    WomSingleFileType.coerceRawValue("/etc/passwd").get shouldEqual WomSingleFile("/etc/passwd")
    WomSingleFileType.coerceRawValue(WomOptionalValue(WomSingleFileType, Option(WomSingleFile("/etc/passwd")))).get shouldEqual WomSingleFile("/etc/passwd")
    WomSingleFileType.coerceRawValue(-1).isFailure shouldBe true
  }

  "WomInteger" should "support expected coercions" in {
    WomIntegerType.coerceRawValue(42).get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue(WomOptionalValue(WomIntegerType, Option(WomInteger(42)))).get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue("42").get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue(JsString("42")).get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue("FAIL").isFailure shouldBe true
  }

  "WomFloatType" should "support expected coercions" in {
    WomFloatType.coerceRawValue(33.3).get shouldEqual WomFloat(33.3)
    WomFloatType.coerceRawValue(WomOptionalValue(WomFloatType, Option(WomFloat(33.3)))).get shouldEqual WomFloat(33.3)
    WomFloatType.coerceRawValue("33.3").get shouldEqual WomFloat(33.3)
    WomFloatType.coerceRawValue(JsString("33.3")).get shouldEqual WomFloat(33.3)
    WomFloatType.coerceRawValue("FAIL").isFailure shouldBe true
  }

  val womValueRawStrings = Table(
    ("WomSource", "WomType"),
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

  behavior of "lowestCommonSubtype"

  // Type A in, Type B in, expected results
  val lcsTestCases: List[(List[WomType], WomType)] = List(
    (List(WomIntegerType, WomStringType), WomStringType),
    (List(WomPairType(WomStringType, WomIntegerType), WomPairType(WomIntegerType, WomStringType)), WomPairType(WomStringType, WomStringType)),
    (List(WomOptionalType(WomIntegerType), WomOptionalType(WomStringType)), WomOptionalType(WomStringType)),
    (List(WomOptionalType(WomIntegerType), WomStringType), WomOptionalType(WomStringType)),
    (List(WomArrayType(WomOptionalType(WomIntegerType)), WomArrayType(WomOptionalType(WomStringType))), WomArrayType(WomOptionalType(WomStringType))),
    (List(WomOptionalType(WomMapType(WomIntegerType, WomStringType)), WomOptionalType(WomMapType(WomStringType, WomIntegerType))), WomOptionalType(WomMapType(WomStringType, WomStringType))),
    (List(
      WomCompositeType(Map(
        "i" -> WomIntegerType,
        "s" -> WomStringType
      )),
      WomCompositeType(Map(
        "a" -> WomStringType,
        "b" -> WomIntegerType
      ))
    ), WomObjectType),
    (List(WomIntegerType, WomFloatType), WomFloatType),
    (List(WomIntegerType, WomBooleanType), WomStringType),
    (List(WomOptionalType(WomMaybeEmptyArrayType(WomSingleFileType)), WomMaybeEmptyArrayType(WomNothingType)), WomOptionalType(WomMaybeEmptyArrayType(WomSingleFileType))),
    (List(WomMaybeEmptyArrayType(WomSingleFileType), WomMaybeEmptyArrayType(WomNothingType)), WomMaybeEmptyArrayType(WomSingleFileType)),
    (List(WomMaybeEmptyArrayType(WomStringType), WomMaybeEmptyArrayType(WomIntegerType), WomMaybeEmptyArrayType(WomNothingType)), WomMaybeEmptyArrayType(WomStringType))
  )

  lcsTestCases foreach { case (types, expectedLcs) =>
    it should s"choose ${expectedLcs.stableName} as the lowest common subtype of [${types.map(_.stableName).mkString(", ")}]" in {
      WomType.lowestCommonSubtype(types) should be(expectedLcs)
      WomType.lowestCommonSubtype(types.reverse) should be(expectedLcs)
      WomType.lowestCommonSubtype(Random.shuffle(types)) should be(expectedLcs)
    }
  }

}
