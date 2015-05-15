package cromwell

import java.io.File

import cromwell.binding.types.{WdlStringType, WdlIntegerType}
import cromwell.binding.values._
import cromwell.binding.{WdlExpressionException, WdlFunctions, WdlExpression}
import cromwell.binding.{WdlExpression, WdlFunctions}
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Success, Failure, Try}

class WdlExpressionSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString

  def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")
  class NoFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = fail("No functions should be called in this test")
  }

  def identifierLookup(String: String): WdlValue = {
    String match {
      case "a" => WdlInteger(1)
      case "b" => WdlInteger(2)
      case "s" => WdlString("s")
      case "o" => WdlObject(Map("key1" -> WdlString("value1"), "key2" -> WdlInteger(9)))
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

  def constEval(exprStr: String): WdlPrimitive = expr(exprStr).evaluate(noLookup, new NoFunctions()).asInstanceOf[Try[WdlPrimitive]].get
  def constEvalError(exprStr: String): Unit = {
    expr(exprStr).evaluate(noLookup, new NoFunctions()).asInstanceOf[Try[WdlPrimitive]] match {
      case Failure(ex) => // Expected
      case Success(value) => fail(s"Operation was supposed to fail, instead I got value: $value")
    }
  }
  def identifierEval(exprStr: String): WdlPrimitive = expr(exprStr).evaluate(identifierLookup, new TestFunctions()).asInstanceOf[Try[WdlPrimitive]].get
  def identifierEvalError(exprStr: String): Unit = {
    expr(exprStr).evaluate(identifierLookup, new TestFunctions()).asInstanceOf[Try[WdlPrimitive]] match {
      case Failure(ex) => // Expected
      case Success(value) => fail(s"Operation was supposed to fail, instead I got value: $value")
    }
  }

  /* Integers */
  "Expression Evaluator with LHS Int" should "Int + Int = Int" in {
    constEval("1+2") shouldEqual WdlInteger(3)
  }
  it should "Int + String = String" in {
    constEval(""" 1 + "String" """) shouldEqual WdlString("1String")
  }
  it should "Int + float = float" in {
    constEval("1+2.3") shouldEqual WdlFloat(3.3.toFloat)
  }
  it should "Int + Boolean = error" in {
    constEvalError("1+true")
  }
  it should "Int - Int = Int" in {
    constEval("3-5") shouldEqual WdlInteger(-2)
  }
  it should "Int - float = float" in {
    constEval("10-6.6") shouldEqual WdlFloat(3.4.toFloat)
  }
  it should "Int - Boolean = error" in {
    constEvalError("1-true")
  }
  it should "Int - String = error" in {
    constEvalError(""" 1-"s"  """)
  }
  it should "Int * Int = Int" in {
    constEval("8 * 7") shouldEqual WdlInteger(56)
  }
  it should "Int * float = float" in {
    constEval("5 * 7.2") shouldEqual WdlFloat(36.toFloat)
  }
  it should "Int * Boolean = error" in {
    constEvalError("1*true")
  }
  it should "Int * String = error" in {
    constEvalError(""" 1*"s"  """)
  }
  it should "Int / Int = Int" in {
    constEval("80 / 6") shouldEqual WdlInteger(13)
  }
  it should "Int / float = float" in {
    constEval("25/2.0") shouldEqual WdlFloat(12.5.toFloat)
  }
  it should "Int / Boolean = error" in {
    constEvalError("1/true")
  }
  it should "Int / String = error" in {
    constEvalError(""" 1/"s"  """)
  }
  it should "Int % Int = Int" in {
    constEval("10 % 3") shouldEqual WdlInteger(1)
  }
  it should "Int & float = float" in {
    constEval("10 % 3.5") shouldEqual WdlFloat(3.0.toFloat)
  }
  it should "Int % Boolean = error" in {
    constEvalError("1%false")
  }
  it should "Int % String = error" in {
    constEvalError(""" 1%"s"  """)
  }
  it should "Int == Int (return true)" in {
    constEval(" 24 == 24 ") shouldEqual WdlBoolean.True
  }
  it should "Int == Int (return false)" in {
    constEval(" 24 == 26 ") shouldEqual WdlBoolean.False
  }
  it should "Int == Boolean (error)" in {
    constEvalError(" 24 == false ")
  }
  it should "Int == String (error)" in {
    constEvalError(""" 1 == "s"  """)
  }
  it should "Int != Int (return true)" in {
    constEval(" 1 != 0 ") shouldEqual WdlBoolean.True
  }
  it should "Int != Int (return false)" in {
    constEval(" 1 != 1 ") shouldEqual WdlBoolean.False
  }
  it should "Int != Boolean (error)" in {
    constEvalError(" 24 != false ")
  }
  it should "Int != String (error)" in {
    constEvalError(""" 1 != "s"  """)
  }
  it should "Int < Int (return false)" in {
    constEval("4 < 3") shouldEqual WdlBoolean.False
  }
  it should "Int < Int (return true)" in {
    constEval("3 < 4") shouldEqual WdlBoolean.True
  }
  it should "Int < float" in {
    constEval("4 < 5.0") shouldEqual WdlBoolean.True
  }
  it should "Int < Boolean (false)" in {
    constEvalError(" 24 < false ")
  }
  it should "Int < String (error)" in {
    constEvalError(""" 1 < "s"  """)
  }
  it should "Int <= Int (return true)" in {
    constEval("3 <= 4") shouldEqual WdlBoolean.True
  }
  it should "Int <= float (return true)" in {
    constEval("3 <= 3.0") shouldEqual WdlBoolean.True
  }
  it should "Int <= Boolean (error)" in {
    constEvalError(" 24 <= false ")
  }
  it should "Int <= String (error)" in {
    constEvalError(""" 1 <= "s"  """)
  }
  it should "Int > Int (return true)" in {
    constEval("4 > 3") shouldEqual WdlBoolean.True
  }
  it should "Int > float (return true)" in {
    constEval("4 > 3.0") shouldEqual WdlBoolean.True
  }
  it should "Int > Boolean (error)" in {
    constEvalError("4 > false")
  }
  it should "Int > String (error)" in {
    constEvalError(""" 1 > "s"  """)
  }
  it should "Int >= Int (return true)" in {
    constEval("4 >= 4") shouldEqual WdlBoolean.True
  }
  it should "Int >= float (return true)" in {
    constEval("4 >= 4.0") shouldEqual WdlBoolean.True
  }
  it should "Int >= Boolean (error)" in {
    constEvalError("4 >= false")
  }
  it should "Int >= String (error)" in {
    constEvalError(""" 1 >= "s"  """)
  }
  it should "-Int + -Int = Int" in {
    constEval("-1 + -3") shouldEqual WdlInteger(-4)
  }

  /* Floats */
  "Expression Evaluator with LHS Float" should "Float + Int = Float" in {
    constEval("1.0+2") shouldEqual WdlFloat(3.toFloat)
  }
  it should "Float + String = String" in {
    constEval(""" 1.0 + "String" """) shouldEqual WdlString("1.0String")
  }
  it should "Float + float = float" in {
    constEval("1.0+2.3") shouldEqual WdlFloat(3.3.toFloat)
  }
  it should "Float + Boolean = error" in {
    constEvalError("1.0+true")
  }
  it should "Float - Int = Float" in {
    constEval("3.0-5") shouldEqual WdlFloat(-2.toFloat)
  }
  it should "Float - Float = float" in {
    constEval("10.0-6.6") shouldEqual WdlFloat(3.4.toFloat)
  }
  it should "Float - Boolean = error" in {
    constEvalError("1.0-true")
  }
  it should "Float - String = error" in {
    constEvalError(""" 1.0-"s"  """)
  }
  it should "Float * Int = Float" in {
    constEval("8.0 * 7") shouldEqual WdlFloat(56.toFloat)
  }
  it should "Float * Float = float" in {
    constEval("5.0 * 7.2") shouldEqual WdlFloat(36.toFloat)
  }
  it should "Float * Boolean = error" in {
    constEvalError("1.0*true")
  }
  it should "Float * String = error" in {
    constEvalError(""" 1.0*"s"  """)
  }
  it should "Float / Int = Float" in {
    constEval("25.0 / 4") shouldEqual WdlFloat(6.25.toFloat)
  }
  it should "Float / float = float" in {
    constEval("25.0/2.0") shouldEqual WdlFloat(12.5.toFloat)
  }
  it should "Float / Boolean = error" in {
    constEvalError("1.0/true")
  }
  it should "Float / 0.0 = error" in {
    constEvalError("1.0/0.0")
  }
  it should "Float / String = error" in {
    constEvalError(""" 1.0/"s"  """)
  }
  it should "Float % Int = Float" in {
    constEval("10.0 % 3") shouldEqual WdlFloat(1.toFloat)
  }
  it should "Float % 0 = error" in {
    constEvalError("10.0 % 0")
  }
  it should "Float & float = float" in {
    constEval("10.0 % 3.5") shouldEqual WdlFloat(3.0.toFloat)
  }
  it should "Float % Boolean = error" in {
    constEvalError("1.0%false")
  }
  it should "Float % String = error" in {
    constEvalError(""" 1.0%"s"  """)
  }
  it should "Float == Int (return true)" in {
    constEval("24.0 == 24 ") shouldEqual WdlBoolean.True
  }
  it should "Float == Float (return true)" in {
    constEval("24.0 == 24.0 ") shouldEqual WdlBoolean.True
  }
  it should "Float == Int (return false)" in {
    constEval("24.0 == 26 ") shouldEqual WdlBoolean.False
  }
  it should "Float == Boolean (error)" in {
    constEvalError("24.0 == false ")
  }
  it should "Float == String (error)" in {
    constEvalError(""" 1.0 == "s"  """)
  }
  it should "Float != Int (return true)" in {
    constEval("1.0 != 0 ") shouldEqual WdlBoolean.True
  }
  it should "Float != Float (return true)" in {
    constEval("1.0 != 0.0 ") shouldEqual WdlBoolean.True
  }
  it should "Float != Int (return false)" in {
    constEval("1.0 != 1 ") shouldEqual WdlBoolean.False
  }
  it should "Float != Boolean (error)" in {
    constEvalError("24.0 != false ")
  }
  it should "Float != String (error)" in {
    constEvalError(""" 1.0 != "s"  """)
  }
  it should "Float < Int (return false)" in {
    constEval("4.0 < 3") shouldEqual WdlBoolean.False
  }
  it should "Float < Float (return false)" in {
    constEval("4.0 < 3.0") shouldEqual WdlBoolean.False
  }
  it should "Float < Int (return true)" in {
    constEval("3.0 < 4") shouldEqual WdlBoolean.True
  }
  it should "Float < Float" in {
    constEval("4.0 < 5.0") shouldEqual WdlBoolean.True
  }
  it should "Float < Boolean (false)" in {
    constEvalError("24.0 < false ")
  }
  it should "Float < String (error)" in {
    constEvalError(""" 1.0 < "s"  """)
  }
  it should "Float <= Int (return true)" in {
    constEval("3.0 <= 4") shouldEqual WdlBoolean.True
  }
  it should "Float <= Float (return true)" in {
    constEval("3.0 <= 3.0") shouldEqual WdlBoolean.True
  }
  it should "Float <= Boolean (error)" in {
    constEvalError("24.0 <= false ")
  }
  it should "Float <= String (error)" in {
    constEvalError(""" 1.0 <= "s"  """)
  }
  it should "Float > Int (return true)" in {
    constEval("4.0 > 3") shouldEqual WdlBoolean.True
  }
  it should "Float > Float (return true)" in {
    constEval("4.0 > 3.0") shouldEqual WdlBoolean.True
  }
  it should "Float > Boolean (error)" in {
    constEvalError("4.0 > false")
  }
  it should "Float > String (error)" in {
    constEvalError(""" 1 > "s"  """)
  }
  it should "Float >= Int (return true)" in {
    constEval("4.0 >= 4") shouldEqual WdlBoolean.True
  }
  it should "Float >= Float (return true)" in {
    constEval("4.0 >= 4.0") shouldEqual WdlBoolean.True
  }
  it should "Float >= Boolean (error)" in {
    constEvalError("4.0 >= false")
  }
  it should "Float >= String (error)" in {
    constEvalError(""" 1.0 >= "s"  """)
  }
  it should "+Float" in {
    constEval("+1.0") shouldEqual WdlFloat(1.0.toFloat)
  }
  it should "-Float" in {
    constEval("-1.0") shouldEqual WdlFloat(-1.0.toFloat)
  }

  /* Booleans */
  "Expression Evaluator with LHS Boolean" should "Boolean + Int = error" in {
    constEvalError("true+2")
  }
  it should "Boolean + String = error" in {
    constEvalError(""" true + "String" """)
  }
  it should "Boolean + float = error" in {
    constEvalError("true+2.3")
  }
  it should "Boolean + Boolean = error" in {
    constEvalError("false+true")
  }
  it should "Boolean - Int = error" in {
    constEvalError("false-5")
  }
  it should "Boolean - Float = error" in {
    constEvalError("false-6.6")
  }
  it should "Boolean - Boolean = error" in {
    constEvalError("true-true")
  }
  it should "Boolean - String = error" in {
    constEvalError(""" true-"s"  """)
  }
  it should "Boolean * Int = error" in {
    constEvalError("false * 7")
  }
  it should "Boolean * Float = error" in {
    constEvalError("false * 7.2")
  }
  it should "Boolean * Boolean = error" in {
    constEvalError("false*true")
  }
  it should "Boolean * String = error" in {
    constEvalError(""" false*"s"  """)
  }
  it should "Boolean / Int = error" in {
    constEvalError("false / 4")
  }
  it should "Boolean / float = error" in {
    constEvalError("false/2.0")
  }
  it should "Boolean / Boolean = error" in {
    constEvalError("false/true")
  }
  it should "Boolean / String = error" in {
    constEvalError(""" true/"s"  """)
  }
  it should "Boolean % Int = error" in {
    constEvalError("true % 3")
  }
  it should "Boolean & float = error" in {
    constEvalError("true % 3.5")
  }
  it should "Boolean % Boolean = error" in {
    constEvalError("false%false")
  }
  it should "Boolean % String = error" in {
    constEvalError(""" true % "s"  """)
  }
  it should "Boolean == Int (error)" in {
    constEvalError("true == 24 ")
  }
  it should "Boolean == Float (error)" in {
    constEvalError("true == 24.0 ")
  }
  it should "Boolean == Boolean (return true)" in {
    constEval("false == false ") shouldEqual WdlBoolean.True
  }
  it should "Boolean == Boolean (return false)" in {
    constEval("false == true ") shouldEqual WdlBoolean.False
  }
  it should "Boolean == String (error)" in {
    constEvalError("""true == "s"  """)
  }
  it should "Boolean != Int (error)" in {
    constEvalError("true != 0 ")
  }
  it should "Boolean != Float (error)" in {
    constEvalError("true != 0.0 ")
  }
  it should "Boolean != Boolean (return true)" in {
    constEval("true != false ") shouldEqual WdlBoolean.True
  }
  it should "Boolean != Boolean (return false)" in {
    constEval("true != true ") shouldEqual WdlBoolean.False
  }
  it should "Boolean != String (error)" in {
    constEvalError("""true != "s"  """)
  }
  it should "Boolean < Int (error)" in {
    constEvalError("true < 3")
  }
  it should "Boolean < Float (error)" in {
    constEvalError("true < 3.0")
  }
  it should "Boolean < Float" in {
    constEvalError("true < 5.0")
  }
  it should "Boolean < Boolean (return false)" in {
    constEval("true < false") shouldEqual WdlBoolean.False
  }
  it should "Boolean < String (error)" in {
    constEvalError("""true < "s"  """)
  }
  it should "Boolean <= Int (error)" in {
    constEvalError("true <= 4")
  }
  it should "Boolean <= Float (error)" in {
    constEvalError("true <= 3.0")
  }
  it should "Boolean <= Boolean (return true)" in {
    constEval("false <= false ") shouldEqual WdlBoolean.True
  }
  it should "Boolean <= Boolean (return false)" in {
    constEval("true <= false") shouldEqual WdlBoolean.False
  }
  it should "Boolean <= String (error)" in {
    constEvalError("""true <= "s"  """)
  }
  it should "Boolean > Int (error)" in {
    constEvalError("true > 3")
  }
  it should "Boolean > Float (error)" in {
    constEvalError("true > 3.0")
  }
  it should "Boolean > Boolean (return true)" in {
    constEval("true > false") shouldEqual WdlBoolean.True
  }
  it should "Boolean > String (error)" in {
    constEvalError(""" 1 > "s"  """)
  }
  it should "Boolean >= Int (error)" in {
    constEvalError("true >= 4")
  }
  it should "Boolean >= Float (error)" in {
    constEvalError("true >= 4.0")
  }
  it should "Boolean >= Boolean (return true)" in {
    constEval("false >= false ") shouldEqual WdlBoolean.True
  }
  it should "Boolean >= Boolean (return false)" in {
    constEval("false >= true ") shouldEqual WdlBoolean.False
  }
  it should "Boolean >= String (error)" in {
    constEvalError("""true >= "s"  """)
  }
  it should "Boolean || Int (error)" in {
    constEvalError("true || 4")
  }
  it should "Boolean || Float (error)" in {
    constEvalError("true || 4.0")
  }
  it should "Boolean || Boolean (return true)" in {
    constEval("false || true ") shouldEqual WdlBoolean.True
  }
  it should "Boolean || Boolean (return false)" in {
    constEval("false || false ") shouldEqual WdlBoolean.False
  }
  it should "Boolean || String (error)" in {
    constEvalError("""true || "s"  """)
  }
  it should "Boolean && Int (error)" in {
    constEvalError("true && 4")
  }
  it should "Boolean && Float (error)" in {
    constEvalError("true && 4.0")
  }
  it should "Boolean && Boolean (return true)" in {
    constEval("false && true ") shouldEqual WdlBoolean.False
  }
  it should "Boolean && Boolean (return false)" in {
    constEval("true && true ") shouldEqual WdlBoolean.True
  }
  it should "Boolean && String (error)" in {
    constEvalError("""true && "s"  """)
  }
  it should "!Boolean (return true)" in {
    constEval("!true") shouldEqual WdlBoolean.False
  }
  it should "!Boolean (return false)" in {
    constEval("!false") shouldEqual WdlBoolean.True
  }

  /* Strings */
  "Expression Evaluator with LHS String" should "String + Int = String" in {
    constEval(""" "String" + 456 """) shouldEqual WdlString("String456")
  }
  it should "String + String = String" in {
    constEval(""" "hello" + " world" """) shouldEqual WdlString("hello world")
  }
  it should "String + Float = String" in {
    constEval(""" "hello" + 2.1 """) shouldEqual WdlString("hello2.1")
  }
  it should "String + Boolean = error" in {
    constEvalError(""" "hello" + true """)
  }
  it should "String == String (return true)" in {
    constEval(""" "hello" == "hello" """) shouldEqual WdlBoolean.True
  }
  it should "String == String (return false)" in {
    constEval(""" "hello" == "hello2" """) shouldEqual WdlBoolean.False
  }
  it should "String == Boolean (error)" in {
    constEvalError(""" "hello" == true """)
  }
  it should "String != String (return false)" in {
    constEval(""" "hello" != "hello" """) shouldEqual WdlBoolean.False
  }
  it should "String != String (return true)" in {
    constEval(""" "hello" != "hello2" """) shouldEqual WdlBoolean.True
  }
  it should "String != Boolean (error)" in {
    constEvalError(""" "hello" != true """)
  }
  it should "String < String (return true)" in {
    constEval(""" "abc" < "def" """) shouldEqual WdlBoolean.True
  }
  it should "String < Boolean (error)" in {
    constEvalError(""" "hello" < true """)
  }
  it should "String <= String (return true)" in {
    constEval(""" "abc" <= "def" """) shouldEqual WdlBoolean.True
  }
  it should "String > String (return false)" in {
    constEval(""" "abc" > "def" """) shouldEqual WdlBoolean.False
  }
  it should "String >= String (return false)" in {
    constEval(""" "abc" >= "def" """) shouldEqual WdlBoolean.False
  }
  it should "String > Boolean (error)" in {
    constEvalError(""" "hello" > true """)
  }

  /* Lookup Function */
  "Expression Evaluator with Look-up function" should "resolve identifier" in {
    identifierEval("a") shouldEqual WdlInteger(1)
  }
  it should "resolve identifier in expression" in {
    identifierEval("a + 10") shouldEqual WdlInteger(11)
  }
  it should "be able to resolve two identifiers as Ints and add them" in {
    identifierEval("a + b") shouldEqual WdlInteger(3)
  }
  it should "be able to resolve one identifer as an Int and another as a String" in {
    identifierEval("s + a") shouldEqual WdlString("s1")
  }

  "Expression Evaluator with File as LHS" should "File + String = File" in {
    WdlFile(new File("/etc")).add(WdlString("/sudoers")).get shouldEqual WdlFile(new File("/etc/sudoers"))
  }
  it should "File + File = File" in {
    WdlFile(new File("/etc")).add(WdlFile(new File("/sudoers"))).get shouldEqual WdlFile(new File("/etc/sudoers"))
  }
  it should "File + Integer (error)" in {
    WdlFile(new File("/etc")).add(WdlInteger(1)) match {
      case Failure(ex) => // Expected
      case _ => fail(s"Operation was supposed to fail, instead I got value: $value")
    }
  }
  it should "File == File (return true)" in {
    WdlFile(new File("/etc")).equals(WdlFile(new File("/etc"))).get shouldEqual WdlBoolean.True
  }
  it should "File == File (return false)" in {
    WdlFile(new File("/etc")).equals(WdlFile(new File("/etc2"))).get shouldEqual WdlBoolean.False
  }
  it should "File == String (return true)" in {
    WdlFile(new File("/etc")).equals(WdlString("/etc")).get shouldEqual WdlBoolean.True
  }
  it should "File == Integer (error)" in {
    WdlFile(new File("/etc")).equals(WdlInteger(1)) match {
      case Failure(ex) => // Expected
      case _ => fail(s"Operation was supposed to fail, instead I got value: $value")
    }
  }
  it should "File != File (return true)" in {
    WdlFile(new File("/etc")).notEquals(WdlFile(new File("/etc2"))).get shouldEqual WdlBoolean.True
  }

  "Expression Evaluator with Object as LHS" should "Lookup object string attribute" in {
    identifierEval("o.key1") shouldEqual WdlString("value1")
  }
  it should "Lookup object integer attribute" in {
    identifierEval("o.key2") shouldEqual WdlInteger(9)
  }
  it should "Error if key doesn't exist" in {
    identifierEvalError("o.key3")
  }

  /* WDL Function Invocations */
  "Expression Evaluator with Function Calls" should "call a function" in {
    identifierEval("b(1)") shouldEqual WdlInteger(2)
  }
  it should "call a function call as part of an expression" in {
    identifierEval("b(1) + 10") shouldEqual WdlInteger(12)
  }
  it should "be able to call a function that appends 2 Strings" in {
    identifierEval(""" append("hello ", "world") """) shouldEqual WdlString("hello world")
  }
  it should "be able to call a function that appends 4 Strings" in {
    identifierEval(""" append("a", "b", "c", "d") """) shouldEqual WdlString("abcd")
  }

  /* Order of Operations */
  "Expression Evaluator with Order-of-Operations" should "prefer multiplication over addition" in {
    constEval("1+2*3") shouldEqual WdlInteger(7)
  }
  it should "prefer addition over equality" in {
    constEval("1+2==3") shouldEqual WdlBoolean.True
  }
  it should "bind parenthesized expressions closer" in {
    constEval("(1+2)*3") shouldEqual WdlInteger(9)
  }
}
