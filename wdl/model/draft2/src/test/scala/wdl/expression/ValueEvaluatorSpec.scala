package wdl.expression

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft2.model
import wdl.draft2.model.WdlExpression
import wdl.draft2.model.expression.{NoFunctions, WdlStandardLibraryFunctions, WdlStandardLibraryFunctionsType}
import wdl.shared.FileSizeLimitationConfig
import wom.OptionalNotSuppliedException
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}

class ValueEvaluatorSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString

  def noLookup(String: String): WomValue = fail("No identifiers should be looked up in this test")

  def noTypeLookup(String: String): WomType = fail("No identifiers should be looked up in this test")

  def identifierLookup(name: String): WomValue = name match {
    case "a" => WomInteger(1)
    case "b" => WomInteger(2)
    case "s" => WomString("s")
    case "array_str" => WomArray(WomArrayType(WomStringType), Seq("foo", "bar", "baz").map(WomString))
    case "map_str_int" => WomMap(WomMapType(WomStringType, WomIntegerType), Map(
      WomString("a") -> WomInteger(0),
      WomString("b") -> WomInteger(1),
      WomString("c") -> WomInteger(2)
    ))
    case "o" => WomObject(Map("key1" -> WomString("value1"), "key2" -> WomInteger(9)))
    case "myPair" => WomPair(WomInteger(3), WomString("hello"))
    case "etc_f" => WomSingleFile("/etc")
    case "etc2_f" => WomSingleFile("/etc2")
    case "etc_s" => WomString("/etc")
    case "sudoers_f" => WomSingleFile("/sudoers")
    case "sudoers_s" => WomString("/sudoers")
    case "f" => WomFloat(0.5F)
    case "t" => WomBoolean(true)
    case "someIntAsString" => WomOptionalValue(WomString("1"))
    case "someFloatAsString" => WomOptionalValue(WomString("0.5"))
    case "someStr" => WomOptionalValue(WomString("someStr"))
    case "someInt" => WomOptionalValue(WomInteger(1))
    case "someBoolean" => WomOptionalValue(WomBoolean(false))
    case "someFloat" => WomOptionalValue(WomFloat(0.5F))
    case "someFile" => WomOptionalValue(WomSingleFile("file"))
    case "noneStr" => WomOptionalValue.none(WomStringType)
    case "noneInt" => WomOptionalValue.none(WomIntegerType)
    case "noneBool" => WomOptionalValue.none(WomBooleanType)
    case "noneFloat" => WomOptionalValue.none(WomFloatType)
    case "noneFile" => WomOptionalValue.none(WomSingleFileType)
  }


  def identifierTypeLookup(name: String): WomType = identifierLookup(name).womType

  class TestValueFunctions extends WdlStandardLibraryFunctions {
    override def globHelper(pattern: String): Seq[String] = throw new UnsupportedOperationException()

    override def readFile(path: String, sizeLimit: Int): String = throw new UnsupportedOperationException()

    override def writeFile(path: String, content: String): Try[WomSingleFile] = Failure(new UnsupportedOperationException())

    override def stdout(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Failure(new UnsupportedOperationException())

    override def stderr(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Failure(new UnsupportedOperationException())

    override def read_json(params: Seq[Try[WomValue]]): Try[WomValue] = Failure(new UnsupportedOperationException())

    override def write_tsv(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Failure(new UnsupportedOperationException())

    override def write_json(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Failure(new UnsupportedOperationException())

    override def size(params: Seq[Try[WomValue]]): Try[WomFloat] = Failure(new UnsupportedOperationException())

    override def length(params: Seq[Try[WomValue]]): Try[WomInteger] = Failure(new UnsupportedOperationException())

    override def sub(params: Seq[Try[WomValue]]): Try[WomString] = Failure(new UnsupportedOperationException())

    override def range(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new UnsupportedOperationException())

    override def transpose(params: Seq[Try[WomValue]]): Try[WomArray] = Failure(new UnsupportedOperationException())

    def b(params: Seq[Try[WomValue]]): Try[WomValue] =
      Success(WomInteger(params.head.asInstanceOf[Try[WomInteger]].get.value + 1))

    def append(params: Seq[Try[WomValue]]): Try[WomValue] =
      Success(WomString(params.map(_.asInstanceOf[Try[WomString]].get.value).mkString))

    override protected def fileSizeLimitationConfig: FileSizeLimitationConfig = throw new UnsupportedOperationException
  }

  class TestTypeFunctions extends WdlStandardLibraryFunctionsType {
    def b(params: Seq[Try[WomType]]): Try[WomType] = Success(WomIntegerType)

    def append(params: Seq[Try[WomType]]): Try[WomType] = Success(WomStringType)
  }

  def constEval(exprStr: String): WomValue = expr(exprStr).evaluate(noLookup, new TestValueFunctions()).get

  def constEvalType(exprStr: String): WomType = expr(exprStr).evaluateType(identifierTypeLookup, new TestTypeFunctions).get

  def constEvalError(exprStr: String): Throwable = {
    expr(exprStr).evaluate(noLookup, new TestValueFunctions()).asInstanceOf[Try[WomPrimitive]] match {
      case Failure(ex) => ex
      case Success(v) => fail(s"Operation was supposed to fail, instead I got value: $v")
    }
  }

  def identifierEval(exprStr: String): WomPrimitive = expr(exprStr).evaluate(identifierLookup, new TestValueFunctions()).asInstanceOf[Try[WomPrimitive]].get

  def identifierEvalError(exprStr: String): Unit = {
    expr(exprStr).evaluate(identifierLookup, new TestValueFunctions()).asInstanceOf[Try[WomPrimitive]] match {
      case Failure(_) => // Expected
      case Success(v) => fail(s"Operation was supposed to fail, instead I got value: $v")
    }
  }

  val constantExpressions = Table(
    ("expression", "value"),

    // Integers
    ("1+2", WomInteger(3)),
    (""" 1 + "String" """, WomString("1String")),
    ("1+2.3", WomFloat(3.3)),
    ("3-5", WomInteger(-2)),
    ("10-6.7", WomFloat(3.3)),
    ("8 * 7", WomInteger(56)),
    ("5 * 7.2", WomFloat(36.toDouble)),
    ("80 / 6", WomInteger(13)),
    ("25/2.0", WomFloat(12.5)),
    ("10 % 3", WomInteger(1)),
    ("10 % 3.5", WomFloat(3.0)),
    (" 24 == 24 ", WomBoolean.True),
    (" 24 == 26 ", WomBoolean.False),
    (" 1 != 0 ", WomBoolean.True),
    (" 1 != 1 ", WomBoolean.False),
    ("4 < 3", WomBoolean.False),
    ("3 < 4", WomBoolean.True),
    ("4 < 5.0", WomBoolean.True),
    ("3 <= 4", WomBoolean.True),
    ("3 <= 3.0", WomBoolean.True),
    ("4 > 3", WomBoolean.True),
    ("4 > 3.0", WomBoolean.True),
    ("4 >= 4", WomBoolean.True),
    ("4 >= 4.0", WomBoolean.True),
    ("-1 + -3", WomInteger(-4)),
    ("+1", WomInteger(1)),

    // Floats
    ("1.0+2", WomFloat(3.toDouble)),
    (""" 1.0 + "String" """, WomString("1.0String")),
    ("1.0+2.3", WomFloat(3.3)),
    ("3.0-5", WomFloat(-2.toDouble)),
    ("10.0-6.7", WomFloat(3.3)),
    ("8.0 * 7", WomFloat(56.toDouble)),
    ("5.0 * 7.2", WomFloat(36.toDouble)),
    ("25.0 / 4", WomFloat(6.25)),
    ("25.0/2.0", WomFloat(12.5)),
    ("10.0 % 3", WomFloat(1.toDouble)),
    ("10.0 % 3.5", WomFloat(3.0)),
    ("24.0 == 24 ", WomBoolean.True),
    ("24.0 == 24.0 ", WomBoolean.True),
    ("24.0 == 26 ", WomBoolean.False),
    ("1.0 != 0 ", WomBoolean.True),
    ("1.0 != 0.0 ", WomBoolean.True),
    ("1.0 != 1 ", WomBoolean.False),
    ("4.0 < 3", WomBoolean.False),
    ("4.0 < 3.0", WomBoolean.False),
    ("3.0 < 4", WomBoolean.True),
    ("4.0 < 5.0", WomBoolean.True),
    ("3.0 <= 4", WomBoolean.True),
    ("3.0 <= 3.0", WomBoolean.True),
    ("4.0 > 3", WomBoolean.True),
    ("4.0 > 3.0", WomBoolean.True),
    ("4.0 >= 4", WomBoolean.True),
    ("4.0 >= 4.0", WomBoolean.True),
    ("+1.0", WomFloat(1.0)),
    ("-1.0", WomFloat(-1.0)),

    // Booleans
    ("false == false ", WomBoolean.True),
    ("false == true ", WomBoolean.False),
    ("true != false ", WomBoolean.True),
    ("true != true ", WomBoolean.False),
    ("true < false", WomBoolean.False),
    ("false <= false ", WomBoolean.True),
    ("true <= false", WomBoolean.False),
    ("true > false", WomBoolean.True),
    ("false >= false ", WomBoolean.True),
    ("false >= true ", WomBoolean.False),
    ("false || true ", WomBoolean.True),
    ("false || false ", WomBoolean.False),
    ("false && true ", WomBoolean.False),
    // Short circuits:
    ("false && 1/0 == 1", WomBoolean.False),
    ("true || 1/0 == 1", WomBoolean.True),
    ("true && true ", WomBoolean.True),
    ("!true", WomBoolean.False),
    ("!false", WomBoolean.True),

    // Strings
    (""" "String" + 456 """, WomString("String456")),
    (""" "hello" + " world" """, WomString("hello world")),
    (""" "hello" + 2.1 """, WomString("hello2.1")),
    (""" "hello" == "hello" """, WomBoolean.True),
    (""" "hello" == "hello2" """, WomBoolean.False),
    (""" "hello" != "hello" """, WomBoolean.False),
    (""" "hello" != "hello2" """, WomBoolean.True),
    (""" "abc" < "def" """, WomBoolean.True),
    (""" "abc" <= "def" """, WomBoolean.True),
    (""" "abc" > "def" """, WomBoolean.False),
    (""" "abc" >= "def" """, WomBoolean.False),

    // Order of Operations
    ("1+2*3", WomInteger(7)),
    ("1+2==3", WomBoolean.True),
    ("(1+2)*3", WomInteger(9)),

    // Simple pair:
    ("(1, \"apple\")", WomPair(WomInteger(1), WomString("apple"))),

    // 1-tuple equivalent to a simple value:
    ("(1)", WomInteger(1)),

    // Array in pair:
    ("(\"hello\", [ 1, 2, 3 ])", WomPair(WomString("hello"), WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(2), WomInteger(3))))),

    // Map to pairs:
    ("""{
    | 1: (1, 2),
    | 2: (2, 3)
    |}
    """.stripMargin, WomMap(WomMapType(WomIntegerType, WomPairType(WomIntegerType, WomIntegerType)), Map(
      WomInteger(1) -> WomPair(WomInteger(1), WomInteger(2)),
      WomInteger(2) -> WomPair(WomInteger(2), WomInteger(3))
    )))
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
    ("false || 4"),
    ("false || 4.0"),
    ("""false || "s"  """),
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
    ("a", WomInteger(1)),
    ("a + 10", WomInteger(11)),
    ("a + b", WomInteger(3)),
    ("s + a", WomString("s1")),
    ("o.key1", WomString("value1")),
    ("o.key2", WomInteger(9)),
    ("myPair.left", WomInteger(3)),
    ("myPair.right", WomString("hello")),

    // Call Functions
    ("b(1)", WomInteger(2)),
    ("b(1) + 10", WomInteger(12)),
    (""" append("hello ", "world") """, WomString("hello world")),
    (""" append("a", "b", "c", "d") """, WomString("abcd")),

    // String Interpolation
    // NB: have to use s"..." interpolation and $$ for these to avoid the "possible missing interpolator" warnings. Sigh...
    (s""""prefix.$${a}.suffix"""", WomString("prefix.1.suffix")),
    (s""""prefix.$${a}$${a}.suffix$${a}"""", WomString("prefix.11.suffix1")),
    (s""""$${b}prefix.$${b}$${a}$${a}.suffix$${a}"""", WomString("2prefix.211.suffix1")),
    (s""""$${s}...$${s}"""", WomString("s...s")),
    (s""""$${someStr}"""", WomString("someStr")),
    (s""""$${noneStr}"""", WomString("")),

    // Array Indexing
    ("array_str[0]", WomString("foo")),
    ("array_str[1]", WomString("bar")),
    ("array_str[2]", WomString("baz")),

    // Map Indexing
    ("""map_str_int["a"]""", WomInteger(0)),
    ("""map_str_int["b"]""", WomInteger(1)),
    ("""map_str_int["c"]""", WomInteger(2)),

    // Files -- 'etc_f', 'etc2_f', and 'sudoers_f' are all variables that resolve to WdlFile
    ("etc_f + sudoers_s", WomSingleFile("/etc/sudoers")),
    (""" "/foo" + etc_f """, WomString("/foo/etc")),
    ("etc_f == etc_f", WomBoolean.True),
    ("etc_f == etc2_f", WomBoolean.False),
    ("etc_f != etc2_f", WomBoolean.True),
    ("etc_f == etc_s", WomBoolean.True),

    // String escaping
    (""" "abc" """, WomString("abc")),
    (""" "a\nb" """, WomString("a\nb")),
    (""" "a\nb\t" """, WomString("a\nb\t")),
    (""" "a\n\"b\t\"" """, WomString("a\n\"b\t\"")),
    (""" "be \u266f or be \u266e, just don't be \u266d" """, WomString("be \u266f or be \u266e, just don't be \u266d")),

    // Optional types
    // String
    ("s + someStr", WomString("ssomeStr")),
    ("s + someInt", WomString("s1")),
    ("s + someFloat", WomString("s0.5")),
    ("s + someFile", WomString("sfile")),
    ("s == someStr", WomBoolean(false)),
    ("s < someStr", WomBoolean(true)),
    ("s > someStr", WomBoolean(false)),
    
    ("someStr + s", WomString("someStrs")),
    ("someInt + s", WomString("1s")),
    ("someFloat + s", WomString("0.5s")),
    ("someFile + s", WomSingleFile("files")),
    ("someStr == s", WomBoolean(false)),
    ("someStr < s", WomBoolean(false)),
    ("someStr > s", WomBoolean(true)),

    // Integer
    ("a + someIntAsString", WomString("11")),
    ("a + someInt", WomInteger(2)),
    ("a * someInt", WomInteger(1)),
    ("a / someInt", WomInteger(1)),
    ("a == someInt", WomBoolean(true)),
    ("a > someInt", WomBoolean(false)),
    ("a < someInt", WomBoolean(false)),

    ("someIntAsString + a", WomString("11")),
    ("someInt + a", WomInteger(2)),
    ("someInt * a", WomInteger(1)),
    ("someInt / a", WomInteger(1)),
    ("someInt == a", WomBoolean(true)),
    ("someInt > a", WomBoolean(false)),
    ("someInt < a", WomBoolean(false)),

    ("-someInt", WomInteger(-1)),
    ("+someInt", WomInteger(1)),

    // Float
    ("f + someFloatAsString", WomString("0.50.5")),
    ("f + someFloat", WomFloat(1)),
    ("f * someFloat", WomFloat(0.25)),
    ("f / someFloat", WomFloat(1)),
    ("f == someFloat", WomBoolean(true)),
    ("f > someFloat", WomBoolean(false)),
    ("f < someFloat", WomBoolean(false)),

    ("someFloatAsString + f", WomString("0.50.5")),
    ("someFloat + f", WomFloat(1)),
    ("someFloat * f", WomFloat(0.25)),
    ("someFloat / f", WomFloat(1)),
    ("someFloat == f", WomBoolean(true)),
    ("someFloat > f", WomBoolean(false)),
    ("someFloat < f", WomBoolean(false)),

    ("-someFloat", WomFloat(-0.5)),
    ("+someFloat", WomFloat(0.5)),

    // Boolean
    ("t == someBoolean", WomBoolean(false)),
    ("t > someBoolean", WomBoolean(true)),
    ("t < someBoolean", WomBoolean(false)),
    ("t && someBoolean", WomBoolean(false)),
    ("t || someBoolean", WomBoolean(true)),

    ("someBoolean == t", WomBoolean(false)),
    ("someBoolean > t", WomBoolean(false)),
    ("someBoolean < t", WomBoolean(true)),
    ("someBoolean && t", WomBoolean(false)),
    ("someBoolean || t", WomBoolean(true)),

    ("!someBoolean", WomBoolean(true)),

    // File
    ("etc_f + someStr", WomSingleFile("/etcsomeStr")),
    ("etc_f == someStr", WomBoolean(false)),
    ("etc_f == someFile", WomBoolean(false)),

    ("someFile == etc_f", WomBoolean(false))
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
    it should s"evaluate $expression into ${value.valueString} (${value.womType.stableName})" in {
      constEval(expression) shouldEqual value
    }
  }

  forAll (constantExpressions) { (expression, value) =>
    it should s"evaluate $expression into type ${value.womType.stableName}" in {
      constEvalType(expression) shouldEqual value.womType
    }
  }

  forAll (identifierLookupExpressions) { (expression, value) =>
    it should s"evaluate $expression into ${value.valueString} (${value.womType.stableName})" in {
      identifierEval(expression) shouldEqual value
    }
  }

  forAll (identifierLookupExpressions) { (expression, value) =>
    it should s"evaluate $expression into type ${value.womType.stableName}" in {
      // need to skip the object expressions because we don't know the types of sub-objects
      if (!expression.startsWith("o.key")) constEvalType(expression) shouldEqual value.womType
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
    WomString("a\nb").toWomString shouldEqual "\"a\\nb\""
    WomString("a\nb\t").toWomString shouldEqual "\"a\\nb\\t\""
    WomString("be \u266f or be \u266e, just don't be \u266d").toWomString shouldEqual "\"be \\u266F or be \\u266E, just don't be \\u266D\""
  }

  "Optional values" should "fail to perform addition with the + operator if the argument is None" in {
    val hello = WomString("hello ")
    val noneWorld = WomOptionalValue.none(WomStringType)
    hello.add(noneWorld) should be(Failure(OptionalNotSuppliedException("+")))
  }

  "Ternary if blocks" should "evaluate only the LHS if the condition is true" in {
    constEval(""" if (5 == 4 + 1) then 6 + 7 else fail() """) should be(WomInteger(13))
  }

  "Ternary if blocks" should "evaluate only the RHS if the condition is false" in {
    constEval(""" if 5 + 6 == 7 then fail() else 14 * 15 """) should be(WomInteger(210))
  }

  "Ternary if blocks" should "fail to evaluate if the condition is not a boolean" in {
    constEvalError(""" if 5 + 6 then 6 + 7 else 14 * 15 """).getMessage should be("'if' expression must be given a boolean argument but got: 11")
  }

  "Ternary if blocks" should "fail to evaluate if the chosen LHS expression fails to evaluate" in {
    constEvalError(""" if (5 == 4 + 1) then fail() else 13 """).getClass.getSimpleName should be("NoSuchMethodException")
  }

  "Ternary if blocks" should "fail to evaluate if the chosen RHS expression fails to evaluate" in {
    constEvalError(""" if (5 == 6 + 1) then 13 else fail() """).getClass.getSimpleName should be("NoSuchMethodException")
  }

  "WdlMaps" should "be coerced to their lowest common WomType" in {
    def lookup(str: String) = str match {
      case "hello" => WomOptionalValue(WomString("bonjour"))
      case _ => fail("Au revoir !")
    }

    val exp = WdlExpression.fromString("""{ "hello": hello, "goodbye": "goodbye" }""")
    val evaluated = exp.evaluate(lookup, NoFunctions)
    evaluated.isSuccess shouldBe true
    evaluated.get.womType shouldBe WomMapType(WomStringType, WomOptionalType(WomStringType))
  }

  it should "evaluate a map literal with values of String? and Int? into a Map[String, String?]" in {
    val str = """ { "i": i_in, "s": s_in } """
    val lookup = Map(
      "i_in" -> WomOptionalValue(WomIntegerType, Some(WomInteger(1))),
      "s_in" -> WomOptionalValue(WomStringType, Some(WomString("two")))
    )
    val exp = model.WdlExpression.fromString(str)

    val expectedMap: WomMap = WomMap(WomMapType(WomStringType, WomOptionalType(WomStringType)), Map (
      WomString("i") -> WomOptionalValue(WomStringType, Some(WomString("1"))),
      WomString("s") -> WomOptionalValue(WomStringType, Some(WomString("two")))
    ))

    val evaluated = exp.evaluate(lookup, NoFunctions)

    evaluated.isSuccess should be(true)
    evaluated.get should be(expectedMap)
  }
}
