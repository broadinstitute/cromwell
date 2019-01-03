package cwl

import common.validation.ErrorOr._
import eu.timepit.refined.refineMV
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.callable.RuntimeEnvironment
import wom.values._
import common.validation.Validation._
import cwl.ExpressionInterpolator.SubstitutionException
import eu.timepit.refined.numeric.Positive
import wom.expression.DefaultSizeIoFunctionSet

class ExpressionInterpolatorSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "ExpressionInterpolator"

  private lazy val expressionLib = Vector(
    better.files.File(getClass.getResource("underscore.js").getPath).contentAsString,
    "var t = function(s) { return _.template(s, {variable: 'data'})({'inputs': inputs}); };"
  )

  private lazy val parameterContext = {
    val runtime = RuntimeEnvironment("out", "tmp", refineMV[Positive](1), 2.0D, 100L, 200L)
    val inputs = Map(
      "cram" -> WomSingleFile("/path/to/my.cram"),
      "file1" -> WomSingleFile("/path/to/my.file.txt")
    )
    ParameterContext(DefaultSizeIoFunctionSet, expressionLib, inputs = inputs, runtimeOption = Option(runtime))
  }

  private def evaluator(string: String): ErrorOr[WomValue] = {
    ExpressionEvaluator.eval(string, parameterContext)
  }

  private val validTests = Table(
    ("description", "expression", "expected"),
    ("an empty string", "", WomString("")),
    ("a space", " ", WomString("")),
    ("two spaces", "  ", WomString("")),
    ("whitespace string", " \n \t", WomString("")),
    ("a regular string", "hello world", WomString("hello world")),
    ("a dollar sign", "$6.28", WomString("$6.28")),
    ("two dollar signs", "$$(6.28)", WomString("$$(6.28)")),
    ("three dollar signs", "$$$(6.28)", WomString("$$6.28")),
    ("a backslash and a dollar sign", """\$(6.28)""", WomString("$(6.28)")),
    ("two backslashes and a dollar sign", """\\$(6.28)""", WomString("""\6.28""")),
    ("a underscore and a dollar sign", "_$(6.28)", WomString("_6.28")),
    ("nested parens", """$(parseInt("6")).$(parseInt(("2")) + "" + parseInt(((("8")))))""", WomString("6.28")),
    ("missing open parens", """$(parseInt("6")).$(parseInt(("2")) + "" + parseInt((("8")))))""", WomString("6.28)")),
    ("quoted open parens", """$("(")""", WomString("(")),
    ("quoted closing parens", """$(")")""", WomString(")")),
    ("prefixed/suffixed nested parens",
      """6.$(parseInt("2"))$(parseInt(("8")) + "" + parseInt(((("3")))))1""",
      WomString("6.2831")),
    ("expressions using libs",
      """$(t("The file is <%= data.inputs.file1.path.split('/').slice(-1)[0] %>\n"))""",
      WomString("The file is my.file.txt\n")),
    ("two expressions", """$(runtime.outdir)/$(inputs.cram.basename)""", WomString("out/my.cram")),
    ("two expressions suffixed", "$(runtime.outdir)/$(inputs.cram.basename).crai", WomString("out/my.cram.crai")),
    ("two expressions prefixed/suffixed",
      "/path/to/$(runtime.outdir)/$(inputs.cram.basename).crai",
      WomString("/path/to/out/my.cram.crai")),
    ("a replacement that must be escaped", """$("$(hello)")""", WomString("$(hello)")),
    ("a string containing a start prefix", """$("$(")""", WomString("$(")),
    ("a replacement returning a json object",
      "$({'output': (inputs.i1 == 'the-default' ? 1 : 2)})",
      WomObject(Map("output" -> WomInteger(2))))
  )

  forAll(validTests) { (description, expression, expected) =>
    it should s"interpolate $description" in {
      ExpressionInterpolator.interpolate(expression, evaluator).toTry.get should be(expected)
    }
  }

  private val substitutionErrorTests = Table(
    ("description", "expression", "expected"),
    ("a missing closing paren", "$(", "0: $("),
    ("a missing quote", """$("))""", """0: $("))"""),
    ("complex missing closing parens",
      """$(parseInt("6")).$(parseInt(("2")) + "" + parseInt(((("8"))))""",
      """1: $(parseInt(("2")) + "" + parseInt(((("8"))))""")
  )

  forAll(substitutionErrorTests) { (description, expression, expected) =>
    it should s"throw a substitution error for $description" in {
      the[SubstitutionException] thrownBy {
        ExpressionInterpolator.interpolate(expression, evaluator)
      } should have message s"Substitution error, unfinished block starting at position $expected"
    }
  }

  /**
    * Tests that sorta work but do not generate the same _exact_ responses as cwltool.
    *
    * See docs at top of [[ExpressionInterpolator]].
    *
    * According to cwltool's expression.py when an interpolated string contains more than just an expression,
    * `json.dumps` should be used to serialize Arrays and Maps, and Maps should sort their keys.
    *
    * We don't even test maps with multiple keys here as WomMap does not guarantee key order.
    */
  private val semiSupportedTests = Table(
    ("description", "expression", "expected"),
    ("an array", """$(["hello", "world"])""", """["hello", "world"]"""),
    ("a map", """$({"hello": "world"})""", """object {hello: "world"}"""),
    ("an array of maps", """$([{"hello": "world"}])""", """[object {hello: "world"}]"""),
    ("a map of arrays", """$({"hello": ["world"]})""", """object {hello: ["world"]}""")
  )

  forAll(semiSupportedTests) { (description, expression, expected) =>
    it should s"semi-interpolate a string plus $description" in {
      val actual = ExpressionInterpolator.interpolate("_" + expression, evaluator).map(_.valueString).toTry.get
      actual should be("_" + expected)
    }
  }
}
