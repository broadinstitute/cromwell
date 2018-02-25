package wom.types

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import spray.json.JsString
import wom.values._

import scala.runtime.ScalaRunTime

class WomTypeSpec extends FlatSpec with Matchers {
  "WomType class" should "stringify WomBoolean to 'Boolean'" in {
    WomBooleanType.toDisplayString shouldEqual "Boolean"
  }
  it should "stringify WomInteger to 'Integer'" in {
    WomIntegerType.toDisplayString shouldEqual "Int"
  }
  it should "stringify WomFloat to 'Float'" in {
    WomFloatType.toDisplayString shouldEqual "Float"
  }
  it should "stringify WomObject to 'Object'" in {
    WomObjectType.toDisplayString shouldEqual "Object"
  }
  it should "stringify WomString to 'String'" in {
    WomStringType.toDisplayString shouldEqual "String"
  }
  it should "stringify WomFile to 'File'" in {
    WomSingleFileType.toDisplayString shouldEqual "File"
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
      "No coercion defined from '0' of type 'Int' to 'Boolean'."
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
      """No coercion defined from '\[0, 1, 2\]' of type 'Array\[Int\?\]' to 'Array\[Int\]\?'."""
    ),
    (
      WomArray(WomArrayType(WomOptionalType(WomIntegerType)), Seq(WomOptionalValue.none(WomIntegerType))),
      WomOptionalType(WomMaybeEmptyArrayType(WomIntegerType)),
      classOf[IllegalArgumentException],
      """No coercion defined from '\[null\]' of type 'Array\[Int\?\]' to 'Array\[Int\]\?'."""
    ),
    (
      WomOptionalValue(WomArray(WomArrayType(WomIntegerType), Seq(
        WomInteger(0),
        WomInteger(1),
        WomInteger(2),
        WomInteger(3),
        WomInteger(4)
      ))),
      WomMaybeEmptyArrayType(WomOptionalType(WomIntegerType)),
      classOf[IllegalArgumentException],
      """No coercion defined from '\[0, 1, 2\]' of type 'Array\[Int\]\?' to 'Array\[Int\?\]'."""
    ),
    (
      WomOptionalValue.none(WomArrayType(WomIntegerType)),
      WomMaybeEmptyArrayType(WomOptionalType(WomIntegerType)),
      classOf[IllegalArgumentException],
      """No coercion defined from 'null' of type 'Array\[Int\]\?' to 'Array\[Int\?\]'."""
    )
  )

  private def describe(any: Any): String = {
    any match {
      case womValue: WomValue => s"wom value ${womValue.toWomString} of type ${womValue.womType.toDisplayString}"
      case _ => s"scala value ${ScalaRunTime.stringOf(any)}"
    }
  }

  forAll(rawValuesCoercedToType) { (rawValue, womType, exceptionClass, exceptionMessage) =>
    it should s"fail coercing ${womType.toDisplayString} from ${describe(rawValue)}" in {
      val exception = womType.coerceRawValue(rawValue).failed.get
      exception.getClass should be(exceptionClass)
      exception.getMessage should fullyMatch regex exceptionMessage
    }
  }

  "WomBoolean" should "support expected coercions" in {
    WomBooleanType.coerceRawValue("true").get shouldEqual WomBoolean.True
    WomBooleanType.coerceRawValue("FALSE").get shouldEqual WomBoolean.False
    WomBooleanType.coerceRawValue(false).get shouldEqual WomBoolean.False

    WomBooleanType.coerceRawValue("I like turtles").isFailure shouldBe true
  }

  "WomString" should "support expected coercions" in {
    WomStringType.coerceRawValue("foo").get shouldEqual WomString("foo")
    WomStringType.coerceRawValue(-1).isFailure shouldBe true
  }

  "WomFile" should "support expected coercions" in {
    WomSingleFileType.coerceRawValue("/etc/passwd").get shouldEqual WomSingleFile("/etc/passwd")
    WomSingleFileType.coerceRawValue(-1).isFailure shouldBe true
  }

  "WomInteger" should "support expected coercions" in {
    WomIntegerType.coerceRawValue(42).get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue("42").get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue(JsString("42")).get shouldEqual WomInteger(42)
    WomIntegerType.coerceRawValue("FAIL").isFailure shouldBe true
  }

  "WomFloatType" should "support expected coercions" in {
    WomFloatType.coerceRawValue(33.3).get shouldEqual WomFloat(33.3)
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
}
