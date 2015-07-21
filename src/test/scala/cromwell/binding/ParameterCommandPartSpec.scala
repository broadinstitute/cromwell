package cromwell.binding

import cromwell.binding.command.ParameterCommandPart
import cromwell.binding.types.{WdlIntegerType, WdlArrayType, WdlStringType}
import cromwell.binding.values.{WdlArray, WdlInteger, WdlString, WdlValue}
import cromwell.parser.BackendType
import cromwell.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}

class ParameterCommandPartSpec extends FlatSpec with Matchers {
  val param1 = ParameterCommandPart(WdlStringType, "name", prefix=None, attributes=Map.empty[String, String])
  val param2 = ParameterCommandPart(WdlStringType, "name", prefix=Some("-p "), attributes=Map.empty[String, String])
  val param3 = ParameterCommandPart(WdlStringType, "name", prefix=Some("-p "), attributes=Map("sep" -> ","), postfixQuantifier=Some("*"))
  val param4 = ParameterCommandPart(WdlIntegerType, "id", prefix=None, attributes=Map("default" -> "1"))

  "command parameter" should "stringify correctly" in {
    param1.toString shouldEqual "${String name}"
  }

  "command parameter instantiation" should "prepend a prefix is specified" in {
    param2.instantiate(Map("name" -> WdlString("foobar"))) shouldEqual "-p foobar"
  }

  it should "combine elements together if * postfix quantifier specified" in {
    val array = WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("foo"), WdlString("bar"), WdlString("baz")))
    param3.instantiate(Map("name" -> array)) shouldEqual "-p foo,bar,baz"
  }

  it should "ignore the default value if a value is specified" in {
    param4.instantiate(Map("id" -> WdlInteger(99))) shouldEqual "99"
  }

  it should "use the default value if a value is not specified" in {
    param4.instantiate(Map.empty[FullyQualifiedName, WdlValue]) shouldEqual "1"
  }

  it should "raise exception if it can't instantiate the parameter" in {
    try {
      param1.instantiate(Map.empty[String, WdlValue])
      fail("Expected an exception")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }

  it should "raise exception if a parameter is a WdlExpression" in {
    try {
      param1.instantiate(Map("name" -> WdlExpression.fromString("1+1")))
      fail("Expected an exception")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }

  it should "raise exception if a parameter not the right type" in {
    try {
      param1.instantiate(Map("name" -> WdlInteger(2)))
      fail("Expected an exception")
    } catch {
      case _: UnsupportedOperationException => // expected
    }
  }

  it should "raise exception if a parameter has a * or + postfix quantifier but no 'sep' attribute set" in {
    try {
      WdlNamespace.load(
        """task test {
          |  command { ./script ${stuff*} }
          |}
        """.stripMargin, BackendType.LOCAL)
      fail("Expected an exception")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "raise exception if a parameter specifies the 'default' attribute but no ? or * postfix quantifier" in {
    try {
      WdlNamespace.load(
        """task test {
          |  command { ./script ${default="x" stuff} }
          |}
        """.stripMargin, BackendType.LOCAL)
      fail("Expected an exception")
    } catch {
      case _: SyntaxError => // expected
    }
  }

  it should "raise exception if a parameter specifies the 'default' attribute but no ? or * postfix quantifier (2)" in {
    try {
      WdlNamespace.load(
        """task test {
          |  command { ./script ${default="x" stuff+} }
          |}
        """.stripMargin, BackendType.LOCAL)
      fail("Expected an exception")
    } catch {
      case _: SyntaxError => // expected
    }
  }
}
