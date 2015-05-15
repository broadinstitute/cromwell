package cromwell

import cromwell.binding.types._
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
}
