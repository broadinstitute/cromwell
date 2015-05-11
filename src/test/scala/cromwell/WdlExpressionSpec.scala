package cromwell

import cromwell.binding.WdlImplicits._
import cromwell.binding.values._
import cromwell.binding.{WdlExpression, WdlFunctions}
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Success, Try}

class WdlExpressionSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString

  def noopLookup(string: String): WdlValue = fail("No identifiers should be looked up in this test")
  class NoopFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = fail("No functions should be called in this test")
  }

  def identifierLookup(string: String): WdlValue = {
    string match {
      case "a" => 1.toWdlValue
      case "b" => 2.toWdlValue
      case "s" => "s".toWdlValue
      case _ => throw new NoSuchElementException(s"$string is not a valid identifier")
    }
  }
  class TestFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = {
      def b(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
        Success(WdlInteger(params.head.asInstanceOf[Try[WdlInteger]].get.value + 1))
      }
      def append(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
        Success(WdlString(params.map(_.asInstanceOf[Try[WdlString]].get.value).mkString))
      }

      name match {
        case "b" => b
        case "append" => append
      }
    }
  }

  def constEval(exprStr: String): WdlPrimitive = expr(exprStr).evaluate(noopLookup, new NoopFunctions()).asInstanceOf[Try[WdlPrimitive]].get
  def identifierEval(exprStr: String): WdlPrimitive = expr(exprStr).evaluate(identifierLookup, new TestFunctions()).asInstanceOf[Try[WdlPrimitive]].get

  "Expression Evaluator" should "be able to add two stinkin' integers" in {
    constEval("1+2") shouldEqual WdlInteger(3)
  }
  it should "be able to add an integer and a string" in {
    constEval(""" 1 + "string" """) shouldEqual WdlString("1string")
  }
  it should "be able to add a string and an integer" in {
    constEval(""" "string" + 456 """) shouldEqual WdlString("string456")
  }
  it should "be able to add two strings" in {
    constEval(""" "hello" + " world" """) shouldEqual WdlString("hello world")
  }
  it should "be able to add two floats" in {
    constEval("1.2 + 3.0") shouldEqual WdlFloat(4.2.toFloat)
  }
  it should "be able to resolve an identifier as an integer" in {
    identifierEval("a + 10") shouldEqual WdlInteger(11)
  }
  it should "be able to resolve two identifiers as integers" in {
    identifierEval("a + b") shouldEqual WdlInteger(3)
  }
  it should "be able to resolve one identifer as an integer and another as a string" in {
    identifierEval("s + a") shouldEqual WdlString("s1")
  }
  it should "be able to call a function" in {
    identifierEval("b(1)") shouldEqual WdlInteger(2)
  }
  it should "be able to call a function and add to the return type" in {
    identifierEval("b(1) + 10") shouldEqual WdlInteger(12)
  }
  it should "be able to call a function that appends 2 strings" in {
    identifierEval(""" append("hello ", "world") """) shouldEqual WdlString("hello world")
  }
  it should "be able to call a function that appends 4 strings" in {
    identifierEval(""" append("a", "b", "c", "d") """) shouldEqual WdlString("abcd")
  }
}
