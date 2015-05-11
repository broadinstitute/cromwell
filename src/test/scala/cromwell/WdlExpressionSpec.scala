package cromwell

import cromwell.binding.types.{WdlStringType, WdlIntegerType}
import cromwell.binding.{WdlValue, WdlFunctions, WdlExpression}
import org.scalatest.{FlatSpec, Matchers}

class WdlExpressionSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString

  def noopLookup(string: String): WdlValue = fail("No identifiers should be looked up in this test")
  class NoopFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = fail("No functions should be called in this test")
  }

  def identifierLookup(string: String): WdlValue = {
    string match {
      case "a" => WdlValue(1, WdlIntegerType)
      case "b" => WdlValue(2, WdlIntegerType)
      case "s" => WdlValue("s", WdlStringType)
      case _ => throw new NoSuchElementException(s"$string is not a valid identifier")
    }
  }
  class TestFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = {
      def b(params: Seq[WdlValue]): WdlValue = {
        new WdlValue(params(0).value.asInstanceOf[Integer] + 1, WdlIntegerType)
      }
      def append(params: Seq[WdlValue]): WdlValue = {
        new WdlValue(params.map(_.value.asInstanceOf[String]).mkString, WdlStringType)
      }

      name match {
        case "b" => b
        case "append" => append
      }
    }
  }

  def constEval(str: String): Any = expr(str).evaluate(noopLookup, new NoopFunctions()).value
  def identifierEval(str: String): Any = expr(str).evaluate(identifierLookup, new TestFunctions()).value

  "Expression Evaluator" should "be able to add two stinkin' integers" in {
    constEval("1+2") shouldEqual 3
  }
  it should "be able to add an integer and a string" in {
    constEval(""" 1 + "string" """) shouldEqual "1string"
  }
  it should "be able to add a string and an integer" in {
    constEval(""" "string" + 456 """) shouldEqual "string456"
  }
  it should "be able to add two strings" in {
    constEval(""" "hello" + " world" """) shouldEqual "hello world"
  }
  it should "be able to add two floats" in {
    constEval("1.2 + 3.0") shouldEqual 4.2.toFloat
  }
  it should "be able to resolve an identifier as an integer" in {
    identifierEval("a + 10") shouldEqual 11
  }
  it should "be able to resolve two identifiers as integers" in {
    identifierEval("a + b") shouldEqual 3
  }
  it should "be able to resolve one identifer as an integer and another as a string" in {
    identifierEval("s + a") shouldEqual "s1"
  }
  it should "be able to call a function" in {
    identifierEval("b(1)") shouldEqual 2
  }
  it should "be able to call a function and add to the return type" in {
    identifierEval("b(1) + 10") shouldEqual 12
  }
  it should "be able to call a function that appends 2 strings" in {
    identifierEval(""" append("hello ", "world") """) shouldEqual "hello world"
  }
  it should "be able to call a function that appends 4 strings" in {
    identifierEval(""" append("a", "b", "c", "d") """) shouldEqual "abcd"
  }
}
