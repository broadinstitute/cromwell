package wdl.draft3.transforms.linking.expression.values

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.ast2wdlom.Ast2WdlomSpec.{fromString, parser}
import wdl.draft3.transforms.ast2wdlom._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.NoIoFunctionSet
import wom.values.WomString

class ValueEvaluatorsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "base value evaluator"

  // In pre-1.1 WDL we use Java-flavored regex

  it should "fail to apply a POSIX-flavor sub expression" in {
    val str = """ sub("aB", "[[:lower:]]", "9") """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedString: WomString = WomString("aB") // no substitution

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedString, Seq.empty)
    }
  }

  it should "evaluate a Java-flavor regex in a sub expression correctly" in {
    val str = """ sub("aB", "\\p{Lower}", "9") """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedString: WomString = WomString("9B")

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedString, Seq.empty)
    }
  }
}
