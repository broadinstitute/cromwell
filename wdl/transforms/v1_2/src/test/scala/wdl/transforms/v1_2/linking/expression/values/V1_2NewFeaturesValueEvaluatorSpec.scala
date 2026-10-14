package wdl.transforms.v1_2.linking.expression.values

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, IdentifierLookup, StringLiteral}
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.v1_2.elements.V1_2ExpressionElement.EnvironmentVariableValue
import wdl.transforms.v1_2.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.v1_2.ast2wdlom._
import wom.expression.NoIoFunctionSet
import wom.types.{WomArrayType, WomCompositeType, WomIntegerType, WomStringType}
import wom.values._

class V1_2NewFeaturesValueEvaluatorSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks {

  behavior of "WDL 1.2 value evaluator for new language features"

  private def evaluate(expression: String, inputs: Map[String, WomValue] = Map.empty) =
    fromString[ExpressionElement](expression, parser.parse_e) match {
      case Right(e) => e.evaluateValue(inputs, NoIoFunctionSet, None).map(_.value)
      case Left(errors) => fail(s"Failed to parse '$expression': ${errors.toList.mkString(", ")}")
    }

  private def ints(values: Int*) = WomArray(WomArrayType(WomIntegerType), values.map(WomInteger.apply))
  private def strings(values: String*) = WomArray(WomArrayType(WomStringType), values.map(WomString.apply))

  private val successfulEvaluations = Table(
    ("expression", "expected value"),
    // Exponentiation
    ("2 ** 10", WomInteger(1024)),
    ("2.0 ** 3", WomFloat(8.0)),
    ("4 ** 0.5", WomFloat(2.0)),
    ("2 * 3 ** 2", WomInteger(18)), // `**` binds tighter than `*`
    ("1 + 2 ** 3 ** 2", WomInteger(65)), // `**` is left-associative: (2 ** 3) ** 2
    // String functions
    (""" find("abc123def456", "[0-9]+") """, WomOptionalValue(WomString("123"))),
    (""" find("abcdef", "[0-9]+") """, WomOptionalValue(WomStringType, None)),
    (""" matches("abc123", "^a.+3$") """, WomBoolean(true)),
    (""" matches("abc123", "^b") """, WomBoolean(false)),
    (""" sub("hello world", "([a-z]+) ([a-z]+)", "\\2 \\1") """, WomString("world hello")),
    (""" sub("cost: 5", "([0-9])", "$\\1") """, WomString("cost: $5")),
    // File functions
    (""" join_paths("/a/b", "c/d.txt") """, WomString("/a/b/c/d.txt")),
    (""" join_paths("/a/b/", ["c", "d.txt"]) """, WomString("/a/b/c/d.txt")),
    (""" join_paths(["a", "b", "c.txt"]) """, WomString("a/b/c.txt")),
    // Array functions
    ("contains([1, 2, 3], 2)", WomBoolean(true)),
    ("""contains(["a", "b"], "c")""", WomBoolean(false)),
    ("contains([1, None], None)", WomBoolean(true)),
    ("chunk([1, 2, 3, 4, 5], 2)",
     WomArray(WomArrayType(WomArrayType(WomIntegerType)), Seq(ints(1, 2), ints(3, 4), ints(5)))
    ),
    ("select_first([None, None], 5)", WomInteger(5)),
    ("select_first([None, 3], 5)", WomInteger(3)),
    // Map, struct and object functions
    ("""contains_key({"a": 1}, "a")""", WomBoolean(true)),
    ("""contains_key({"a": 1}, "b")""", WomBoolean(false)),
    ("""contains_key(object { a: 1 }, "a")""", WomBoolean(true)),
    ("""contains_key({"a": {"b": 1}}, ["a", "b"])""", WomBoolean(true)),
    ("""contains_key({"a": {"b": 1}}, ["a", "c"])""", WomBoolean(false)),
    ("""values({"a": 1, "b": 2})""", ints(1, 2)),
    ("""keys(object { a: 1 })""", strings("a")),
    ("""length({"a": 1, "b": 2})""", WomInteger(2)),
    ("""length("hello")""", WomInteger(5)),
    ("""length([1, 2, 3])""", WomInteger(3)),
    // Literals and operators
    (".5 + 1", WomFloat(1.5)),
    ("1e3", WomFloat(1000.0)),
    ("[1, 2, 3][3 - 1]", WomInteger(3)),
    ("\"~{1 == 1}\"", WomString("true")),
    ("\"~{3.5}\"", WomString("3.5")),
    ("select_first([], 5)", WomInteger(5)),
    // Equality with None
    ("None == None", WomBoolean(true)),
    ("None != None", WomBoolean(false)),
    ("select_first([1]) == None", WomBoolean(false))
  )

  forAll(successfulEvaluations) { (expression, expected) =>
    it should s"evaluate '${expression.trim}'" in {
      evaluate(expression) shouldBeValid expected
    }
  }

  private val failedEvaluations = Table(
    ("expression", "expected error"),
    ("2 ** -1", "an Int cannot be raised to a negative Int power"),
    ("chunk([1, 2], 0)", "chunk size must be greater than 0"),
    (""" join_paths("/a", "/b") """, "only the first path may be absolute")
  )

  forAll(failedEvaluations) { (expression, expectedError) =>
    it should s"fail to evaluate '${expression.trim}'" in {
      evaluate(expression) match {
        case cats.data.Validated.Invalid(errors) => errors.toList.mkString should include(expectedError)
        case valid => fail(s"Expected an evaluation failure but got $valid")
      }
    }
  }

  it should "return the keys of a struct in definition order" in {
    val structType = WomCompositeType(Map("first" -> WomStringType, "last" -> WomStringType), Option("Name"))
    val name = WomObject(Map("last" -> WomString("Doe"), "first" -> WomString("Jane")), structType)
    evaluate("keys(name)", Map("name" -> name)) shouldBeValid strings(structType.typeMap.keys.toSeq: _*)
  }

  // Examples from https://github.com/openwdl/wdl/blob/v1.2.1/SPEC.md#multi-line-strings
  private val multilineStrings = Table(
    ("multi-line string", "expected value"),
    ("<<<hello  world>>>", "hello  world"),
    ("<<<   hello  world   >>>", "hello  world"),
    ("<<<   \n    hello  world>>>", "hello  world"),
    ("<<<   \n    hello  world\n    >>>", "hello  world"),
    ("<<<   \n    hello  world\n>>>", "hello  world"),
    ("<<<\n    hello  \\\n        world\n>>>", "hello  world"),
    ("<<<\nhello \\\\\n  world\n>>>", "hello \\\n  world"),
    ("<<<\n\n      this is a\n      \n        multi-line string\n      \n      >>>",
     "\nthis is a\n\n  multi-line string\n"
    ),
    ("<<<\n      \n        this is a\n        \n          multi-line string\n\n      >>>",
     "\nthis is a\n\n  multi-line string\n"
    ),
    ("<<<\n        multi-line string \\\n        with 'single' and \"double\" quotes\n      >>>",
     "multi-line string with 'single' and \"double\" quotes"
    ),
    ("<<<\n    tab\\tescape\n    >>>", "tab\tescape")
  )

  forAll(multilineStrings) { (multiline, expected) =>
    it should s"evaluate the multi-line string ${multiline.replace("\n", "\\n")}" in {
      evaluate(multiline) shouldBeValid WomString(expected)
    }
  }

  it should "evaluate placeholders in a multi-line string" in {
    val multiline = "<<<\n    hello ~{name}\n      and ~{1 + 1}\n  >>>"
    evaluate(multiline, Map("name" -> WomString("world"))) shouldBeValid WomString("hello world\n  and 2")
  }

  it should "compare optional values with None" in {
    val inputs = Map(
      "i" -> WomInteger(1),
      "j" -> WomOptionalValue(WomInteger(1)),
      "k" -> WomOptionalValue(WomIntegerType, None)
    )
    evaluate("i == j", inputs) shouldBeValid WomBoolean(true)
    evaluate("k == None", inputs) shouldBeValid WomBoolean(true)
    evaluate("i == k", inputs) shouldBeValid WomBoolean(false)
    evaluate("j == k", inputs) shouldBeValid WomBoolean(false)
    evaluate("j != k", inputs) shouldBeValid WomBoolean(true)
  }

  it should "single-quote environment variable values so the shell does not interpret them" in {
    val element: ExpressionElement = EnvironmentVariableValue(StringLiteral("it's $HOME; rm -rf /"))
    element.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(
      WomString("""'it'"'"'s $HOME; rm -rf /'"""),
      Seq.empty
    )
  }

  it should "serialize compound environment variable values as JSON" in {
    val element: ExpressionElement =
      EnvironmentVariableValue(ArrayLiteral(Seq(StringLiteral("a"), StringLiteral("b"))))
    element.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomString("""'["a","b"]'"""),
                                                                                         Seq.empty
    )
  }

  it should "export None environment variable values as an empty string" in {
    val element: ExpressionElement = EnvironmentVariableValue(IdentifierLookup("maybe"))
    element.evaluateValue(Map("maybe" -> WomOptionalValue(WomStringType, None)), NoIoFunctionSet, None) shouldBeValid
      EvaluatedValue(WomString("''"), Seq.empty)
  }
}
