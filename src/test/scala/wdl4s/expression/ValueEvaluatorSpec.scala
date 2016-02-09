package wdl4s.expression

import wdl4s.WdlExpression
import wdl4s.types._
import wdl4s.values._
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Failure, Success, Try}

class ValueEvaluatorSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString

  def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")
  def noTypeLookup(String: String): WdlType = fail("No identifiers should be looked up in this test")

  def identifierLookup(name: String): WdlValue = name match {
    case "a" => WdlInteger(1)
    case "b" => WdlInteger(2)
    case "s" => WdlString("s")
    case "array_str" => WdlArray(WdlArrayType(WdlStringType), Seq("foo", "bar", "baz").map(WdlString))
    case "map_str_int" => WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
      WdlString("a") -> WdlInteger(0),
      WdlString("b") -> WdlInteger(1),
      WdlString("c") -> WdlInteger(2)
    ))
    case "o" => WdlObject(Map("key1" -> WdlString("value1"), "key2" -> WdlInteger(9)))
    case "etc_f" => WdlFile("/etc")
    case "etc2_f" => WdlFile("/etc2")
    case "etc_s" => WdlString("/etc")
    case "sudoers_f" => WdlFile("/sudoers")
    case "sudoers_s" => WdlString("/sudoers")
  }


  def identifierTypeLookup(name: String): WdlType = identifierLookup(name).wdlType

  class TestValueFunctions extends WdlStandardLibraryFunctions {
    override def glob(path: String, pattern: String): Seq[String] = throw new NotImplementedError()
    override def readFile(path: String): String = throw new NotImplementedError()
    override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = throw new NotImplementedError()
    override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
    override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
    override def read_json(params: Seq[Try[WdlValue]]): Try[WdlValue] = Failure(new NotImplementedError())
    override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
    override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError())
    override def size(params: Seq[Try[WdlValue]]): Try[WdlFloat] = Failure(new NotImplementedError())

    def b(params: Seq[Try[WdlValue]]): Try[WdlValue] =
      Success(WdlInteger(params.head.asInstanceOf[Try[WdlInteger]].get.value + 1))

    def append(params: Seq[Try[WdlValue]]): Try[WdlValue] =
      Success(WdlString(params.map(_.asInstanceOf[Try[WdlString]].get.value).mkString))
  }

  class TestTypeFunctions extends WdlStandardLibraryFunctionsType {
    def b(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlIntegerType)
    def append(params: Seq[Try[WdlType]]): Try[WdlType] = Success(WdlStringType)
  }

  def constEval(exprStr: String): WdlPrimitive = expr(exprStr).evaluate(noLookup, new TestValueFunctions()).asInstanceOf[Try[WdlPrimitive]].get
  def constEvalType(exprStr: String): WdlType = expr(exprStr).evaluateType(identifierTypeLookup, new TestTypeFunctions).asInstanceOf[Try[WdlType]].get
  def constEvalError(exprStr: String): Unit = {
    expr(exprStr).evaluate(noLookup, new TestValueFunctions()).asInstanceOf[Try[WdlPrimitive]] match {
      case Failure(ex) => // Expected
      case Success(value) => fail(s"Operation was supposed to fail, instead I got value: $value")
    }
  }
  def identifierEval(exprStr: String): WdlPrimitive = expr(exprStr).evaluate(identifierLookup, new TestValueFunctions()).asInstanceOf[Try[WdlPrimitive]].get
  def identifierEvalError(exprStr: String): Unit = {
    expr(exprStr).evaluate(identifierLookup, new TestValueFunctions()).asInstanceOf[Try[WdlPrimitive]] match {
      case Failure(ex) => // Expected
      case Success(value) => fail(s"Operation was supposed to fail, instead I got value: $value")
    }
  }

  val constantExpressions = Table(
    ("expression", "value"),

    // Integers
    ("1+2", WdlInteger(3)),
    (""" 1 + "String" """, WdlString("1String")),
    ("1+2.3", WdlFloat(3.3)),
    ("3-5", WdlInteger(-2)),
    ("10-6.7", WdlFloat(3.3)),
    ("8 * 7", WdlInteger(56)),
    ("5 * 7.2", WdlFloat(36.toDouble)),
    ("80 / 6", WdlInteger(13)),
    ("25/2.0", WdlFloat(12.5)),
    ("10 % 3", WdlInteger(1)),
    ("10 % 3.5", WdlFloat(3.0)),
    (" 24 == 24 ", WdlBoolean.True),
    (" 24 == 26 ", WdlBoolean.False),
    (" 1 != 0 ", WdlBoolean.True),
    (" 1 != 1 ", WdlBoolean.False),
    ("4 < 3", WdlBoolean.False),
    ("3 < 4", WdlBoolean.True),
    ("4 < 5.0", WdlBoolean.True),
    ("3 <= 4", WdlBoolean.True),
    ("3 <= 3.0", WdlBoolean.True),
    ("4 > 3", WdlBoolean.True),
    ("4 > 3.0", WdlBoolean.True),
    ("4 >= 4", WdlBoolean.True),
    ("4 >= 4.0", WdlBoolean.True),
    ("-1 + -3", WdlInteger(-4)),
    ("+1", WdlInteger(1)),

    // Floats
    ("1.0+2", WdlFloat(3.toDouble)),
    (""" 1.0 + "String" """, WdlString("1.0String")),
    ("1.0+2.3", WdlFloat(3.3)),
    ("3.0-5", WdlFloat(-2.toDouble)),
    ("10.0-6.7", WdlFloat(3.3)),
    ("8.0 * 7", WdlFloat(56.toDouble)),
    ("5.0 * 7.2", WdlFloat(36.toDouble)),
    ("25.0 / 4", WdlFloat(6.25)),
    ("25.0/2.0", WdlFloat(12.5)),
    ("10.0 % 3", WdlFloat(1.toDouble)),
    ("10.0 % 3.5", WdlFloat(3.0)),
    ("24.0 == 24 ", WdlBoolean.True),
    ("24.0 == 24.0 ", WdlBoolean.True),
    ("24.0 == 26 ", WdlBoolean.False),
    ("1.0 != 0 ", WdlBoolean.True),
    ("1.0 != 0.0 ", WdlBoolean.True),
    ("1.0 != 1 ", WdlBoolean.False),
    ("4.0 < 3", WdlBoolean.False),
    ("4.0 < 3.0", WdlBoolean.False),
    ("3.0 < 4", WdlBoolean.True),
    ("4.0 < 5.0", WdlBoolean.True),
    ("3.0 <= 4", WdlBoolean.True),
    ("3.0 <= 3.0", WdlBoolean.True),
    ("4.0 > 3", WdlBoolean.True),
    ("4.0 > 3.0", WdlBoolean.True),
    ("4.0 >= 4", WdlBoolean.True),
    ("4.0 >= 4.0", WdlBoolean.True),
    ("+1.0", WdlFloat(1.0)),
    ("-1.0", WdlFloat(-1.0)),

    // Booleans
    ("false == false ", WdlBoolean.True),
    ("false == true ", WdlBoolean.False),
    ("true != false ", WdlBoolean.True),
    ("true != true ", WdlBoolean.False),
    ("true < false", WdlBoolean.False),
    ("false <= false ", WdlBoolean.True),
    ("true <= false", WdlBoolean.False),
    ("true > false", WdlBoolean.True),
    ("false >= false ", WdlBoolean.True),
    ("false >= true ", WdlBoolean.False),
    ("false || true ", WdlBoolean.True),
    ("false || false ", WdlBoolean.False),
    ("false && true ", WdlBoolean.False),
    ("true && true ", WdlBoolean.True),
    ("!true", WdlBoolean.False),
    ("!false", WdlBoolean.True),

    // Strings
    (""" "String" + 456 """, WdlString("String456")),
    (""" "hello" + " world" """, WdlString("hello world")),
    (""" "hello" + 2.1 """, WdlString("hello2.1")),
    (""" "hello" == "hello" """, WdlBoolean.True),
    (""" "hello" == "hello2" """, WdlBoolean.False),
    (""" "hello" != "hello" """, WdlBoolean.False),
    (""" "hello" != "hello2" """, WdlBoolean.True),
    (""" "abc" < "def" """, WdlBoolean.True),
    (""" "abc" <= "def" """, WdlBoolean.True),
    (""" "abc" > "def" """, WdlBoolean.False),
    (""" "abc" >= "def" """, WdlBoolean.False),

    // Order of Operations
    ("1+2*3", WdlInteger(7)),
    ("1+2==3", WdlBoolean.True),
    ("(1+2)*3", WdlInteger(9))
  )

  val badExpressions = Table(
    ("expression"),

    // Integers
    ("1+true"),
    ("1-true"),
    (""" 1-"s"  """),
    ("1*true"),
    (""" 1*"s"  """),
    ("1 / 0"),
    ("1 / 0.0"),
    ("25/0.0"),
    ("1/true"),
    (""" 1/"s"  """),
    ("1%false"),
    (""" 1%"s"  """),
    ("1%0"),
    (" 24 == false "),
    (""" 1 == "s"  """),
    (" 24 != false "),
    (""" 1 != "s"  """),
    (" 24 < false "),
    (""" 1 < "s"  """),
    (" 24 <= false "),
    (""" 1 <= "s"  """),
    ("4 > false"),
    (""" 1 > "s"  """),
    ("4 >= false"),
    (""" 1 >= "s"  """),

    // Floats
    ("1.0+true"),
    ("1.0-true"),
    (""" 1.0-"s"  """),
    ("1.0*true"),
    (""" 1.0*"s"  """),
    ("1.0/true"),
    ("1.0/0.0"),
    ("1.0/0"),
    (""" 1.0/"s"  """),
    ("10.0 % 0"),
    ("10.0 % 0.0"),
    ("1.0%false"),
    (""" 1.0%"s"  """),
    ("24.0 == false "),
    (""" 1.0 == "s"  """),
    ("24.0 != false "),
    (""" 1.0 != "s"  """),
    ("24.0 < false "),
    (""" 1.0 < "s"  """),
    ("24.0 <= false "),
    (""" 1.0 <= "s"  """),
    ("4.0 > false"),
    (""" 1.0 > "s"  """),
    ("4.0 >= false"),
    (""" 1.0 >= "s"  """),

    // Booleans
    (""" true + "String" """),
    ("true+2"),
    ("true+2.3"),
    ("false+true"),
    ("false-5"),
    ("false-6.6"),
    ("true-true"),
    (""" true-"s"  """),
    ("false * 7"),
    ("false * 7.2"),
    ("false*true"),
    (""" false*"s"  """),
    ("false / 4"),
    ("false/2.0"),
    ("false/true"),
    (""" true/"s"  """),
    ("true % 3"),
    ("true % 3.5"),
    ("false%false"),
    (""" true % "s"  """),
    ("true == 24 "),
    ("true == 24.0 "),
    ("""true == "s"  """),
    ("true != 0 "),
    ("true != 0.0 "),
    ("""true != "s"  """),
    ("true < 3"),
    ("true < 3.0"),
    ("true < 5.0"),
    ("""true < "s"  """),
    ("true <= 4"),
    ("true <= 3.0"),
    ("""true <= "s"  """),
    ("true > 3"),
    ("true > 3.0"),
    ("true >= 4"),
    ("true >= 4.0"),
    ("""true >= "s"  """),
    ("true || 4"),
    ("true || 4.0"),
    ("""true || "s"  """),
    ("true && 4"),
    ("true && 4.0"),
    ("""true && "s"  """),

    // Strings
    (""" "hello" + true """),
    (""" "hello" == true """),
    (""" "hello" != true """),
    (""" "hello" < true """),
    (""" "hello" > true """)
  )

  val identifierLookupExpressions = Table(
    ("expression", "value"),

    // Lookup Variables
    ("a", WdlInteger(1)),
    ("a + 10", WdlInteger(11)),
    ("a + b", WdlInteger(3)),
    ("s + a", WdlString("s1")),
    ("o.key1", WdlString("value1")),
    ("o.key2", WdlInteger(9)),

    // Call Functions
    ("b(1)", WdlInteger(2)),
    ("b(1) + 10", WdlInteger(12)),
    (""" append("hello ", "world") """, WdlString("hello world")),
    (""" append("a", "b", "c", "d") """, WdlString("abcd")),

    // String Interpolation
    ("\"prefix.${a}.suffix\"", WdlString("prefix.1.suffix")),
    ("\"prefix.${a}${a}.suffix${a}\"", WdlString("prefix.11.suffix1")),
    ("\"${b}prefix.${b}${a}${a}.suffix${a}\"", WdlString("2prefix.211.suffix1")),
    ("\"${s}...${s}\"", WdlString("s...s")),

    // Array Indexing
    ("array_str[0]", WdlString("foo")),
    ("array_str[1]", WdlString("bar")),
    ("array_str[2]", WdlString("baz")),

    // Map Indexing
    ("""map_str_int["a"]""", WdlInteger(0)),
    ("""map_str_int["b"]""", WdlInteger(1)),
    ("""map_str_int["c"]""", WdlInteger(2)),

    // Files -- 'etc_f', 'etc2_f', and 'sudoers_f' are all variables that resolve to WdlFile
    ("etc_f + sudoers_s", WdlFile("/etc/sudoers")),
    (""" "/foo" + etc_f """, WdlString("/foo/etc")),
    ("etc_f == etc_f", WdlBoolean.True),
    ("etc_f == etc2_f", WdlBoolean.False),
    ("etc_f != etc2_f", WdlBoolean.True),
    ("etc_f == etc_s", WdlBoolean.True),

    // String escaping
    (""" "abc" """, WdlString("abc")),
    (""" "a\nb" """, WdlString("a\nb")),
    (""" "a\nb\t" """, WdlString("a\nb\t")),
    (""" "a\n\"b\t\"" """, WdlString("a\n\"b\t\"")),
    (""" "be \u266f or be \u266e, just don't be \u266d" """, WdlString("be \u266f or be \u266e, just don't be \u266d"))
  )

  val badIdentifierExpressions = Table(
    ("expression"),
    ("etc_f + 1"),
    ("etc_f == 1"),
    ("0.key3"),
    ("array_str[3]"),
    ("""map_str_int["d"]""")
  )

  forAll (constantExpressions) { (expression, value) =>
    it should s"evaluate $expression into ${value.valueString} (${value.wdlType.toWdlString})" in {
      constEval(expression) shouldEqual value
    }
  }

  forAll (constantExpressions) { (expression, value) =>
    it should s"evaluate $expression into type ${value.wdlType.toWdlString}" in {
      constEvalType(expression) shouldEqual value.wdlType
    }
  }

  forAll (identifierLookupExpressions) { (expression, value) =>
    it should s"evaluate $expression into ${value.valueString} (${value.wdlType.toWdlString})" in {
      identifierEval(expression) shouldEqual value
    }
  }

  forAll (identifierLookupExpressions) { (expression, value) =>
    it should s"evaluate $expression into type ${value.wdlType.toWdlString}" in {
      // need to skip the object expressions because we don't know the types of sub-objects
      if (!expression.startsWith("o.key")) constEvalType(expression) shouldEqual value.wdlType
    }
  }

  forAll (badExpressions) { (expression) =>
    it should s"error when evaluating: $expression" in {
      constEvalError(expression)
    }
  }

  forAll (badIdentifierExpressions) { (expression) =>
    it should s"error when evaluating: $expression" in {
      identifierEvalError(expression)
    }
  }

  "A string with special characters in it" should "convert to escape sequences when converted to WDL" in {
    WdlString("a\nb").toWdlString shouldEqual "\"a\\nb\""
    WdlString("a\nb\t").toWdlString shouldEqual "\"a\\nb\\t\""
    WdlString("be \u266f or be \u266e, just don't be \u266d").toWdlString shouldEqual "\"be \\u266F or be \\u266E, just don't be \\u266D\""
  }
}
