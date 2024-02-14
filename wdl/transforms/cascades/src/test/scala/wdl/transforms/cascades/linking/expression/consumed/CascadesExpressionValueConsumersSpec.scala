package wdl.transforms.cascades.linking.expression.consumed

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedIdentifierHook}
import wdl.transforms.cascades.Ast2WdlomSpec._
import wdl.transforms.cascades.ast2wdlom._

class CascadesExpressionValueConsumersSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  behavior of "WDL Cascades' expressionValueConsumer"

  it should "return nothing from static integer addition" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set.empty)
    }
  }

  it should "find the variable lookup within an integer addition" in {
    val str = "3 + three"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("three")))
    }
  }

  it should "find the variable lookup within as_map(...) called on a task output" in {
    val str = "as_map(my_task.out)"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedCallOutputOrIdentifierAndMemberAccessHook("my_task", "out")))
    }
  }

  it should "find the variable lookup within as_pairs(as_map(...)), called on a task output" in {
    val str = "as_pairs(as_map(my_task.out))"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedCallOutputOrIdentifierAndMemberAccessHook("my_task", "out")))
    }
  }

  it should "discover the variable lookups within a sep() call" in {
    val str = """ sep(my_separator, ["a", "b", c]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_separator"), UnlinkedIdentifierHook("c")))
    }
  }

  it should "discover the variable lookups within a suffix() call" in {
    val str = """ suffix(my_suffix, ["a", "b", c]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_suffix"), UnlinkedIdentifierHook("c")))
    }
  }

  it should "discover an array variable lookup within a suffix() call" in {
    val str = """ suffix("SFX", my_array) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_array")))
    }
  }
}
