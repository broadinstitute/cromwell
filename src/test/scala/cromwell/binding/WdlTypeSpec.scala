package cromwell.binding

import java.nio.file.Paths

import cromwell.binding.types._
import cromwell.binding.values._
import org.scalatest.{FlatSpec, Matchers}

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
    WdlIntegerType.coerceRawValue("FAIL").isFailure shouldBe true
  }

  "WdlFloatType" should "support expected coercions" in {
    WdlFloatType.coerceRawValue(33.3).get shouldEqual WdlFloat(33.3)
    WdlFloatType.coerceRawValue("FAIL").isFailure shouldBe true
  }
}
