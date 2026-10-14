package wdl.transforms.v1_2.linking.expression.types

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.transforms.v1_2.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.v1_2.ast2wdlom._
import wom.types._

class V1_2NewFeaturesTypeEvaluatorSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks {

  behavior of "WDL 1.2 type evaluator for new language features"

  private def evaluateType(expression: String) =
    fromString[ExpressionElement](expression, parser.parse_e) match {
      case Right(e) => e.evaluateType(Map.empty, Map.empty)
      case Left(errors) => fail(s"Failed to parse '$expression': ${errors.toList.mkString(", ")}")
    }

  private val validExpressions = Table(
    ("expression", "expected type"),
    ("2 ** 3", WomIntegerType),
    ("2 ** 0.5", WomFloatType),
    ("2.0 ** 2", WomFloatType),
    (""" find("abc", "b") """, WomOptionalType(WomStringType)),
    (""" matches("abc", "b") """, WomBooleanType),
    (""" join_paths("/a", "b") """, WomStringType),
    (""" join_paths("/a", ["b", "c"]) """, WomStringType),
    (""" join_paths(["/a", "b"]) """, WomStringType),
    (""" read_tsv("f.tsv", true) """, WomArrayType(WomObjectType)),
    (""" read_tsv("f.tsv", false) """, WomArrayType(WomArrayType(WomStringType))),
    (""" read_tsv("f.tsv", false, ["a", "b"]) """, WomArrayType(WomObjectType)),
    (""" write_tsv([["a", "b"]], true, ["x", "y"]) """, WomSingleFileType),
    ("contains([1, 2], 1)", WomBooleanType),
    ("chunk([1, 2, 3], 2)", WomArrayType(WomArrayType(WomIntegerType))),
    ("select_first([1, None], 5)", WomIntegerType),
    ("select_first([], 5)", WomIntegerType),
    ("""contains_key({"a": 1}, "a")""", WomBooleanType),
    ("""contains_key(object { a: 1 }, ["a", "b"])""", WomBooleanType),
    ("""values({"a": 1})""", WomArrayType(WomIntegerType)),
    ("""keys(object { a: 1 })""", WomArrayType(WomStringType)),
    ("""length({"a": 1})""", WomIntegerType),
    ("""length("abc")""", WomIntegerType),
    ("None == None", WomBooleanType),
    ("select_first([1]) == None", WomBooleanType),
    (".5 + 1", WomFloatType),
    ("1e3", WomFloatType),
    ("[1, 2, 3][3 - 1]", WomIntegerType),
    // `read_object` and friends are deprecated in WDL 1.2 but available again
    (""" read_object("f.tsv") """, WomObjectType)
  )

  forAll(validExpressions) { (expression, expectedType) =>
    it should s"evaluate the type of '${expression.trim}'" in {
      evaluateType(expression) shouldBeValid expectedType
    }
  }

  private val invalidExpressions = Table(
    ("expression", "expected error"),
    (""" "a" ** 2 """, "Both operands must be numeric"),
    ("""contains([1, 2], [1])""", "Cannot call 'contains'"),
    ("""values([1, 2])""", "Expected a map"),
    ("""length(1)""", "Expected an array, map, object or string"),
    (""" join_paths({"a": 1}, "a") """, "Invalid first argument to 'join_paths'")
  )

  forAll(invalidExpressions) { (expression, expectedError) =>
    it should s"fail to evaluate the type of '${expression.trim}'" in {
      evaluateType(expression) match {
        case cats.data.Validated.Invalid(errors) => errors.toList.mkString should include(expectedError)
        case valid => fail(s"Expected a type evaluation failure but got $valid")
      }
    }
  }
}
