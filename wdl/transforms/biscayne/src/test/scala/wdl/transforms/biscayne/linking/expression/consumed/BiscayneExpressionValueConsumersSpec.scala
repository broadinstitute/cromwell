package wdl.transforms.biscayne.linking.expression.consumed

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.{UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedIdentifierHook}
import wdl.transforms.biscayne.Ast2WdlomSpec._
import wdl.transforms.biscayne.ast2wdlom._

class BiscayneExpressionValueConsumersSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  behavior of "WDL Biscayne's expressionValueConsumer"

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

  it should "discover the variable lookups within a sub() call" in {
    val str = """ sub(my_input, "^[A-Z]$", "0") """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_input")))
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

  it should "discover an array variable lookup within a quote() call" in {
    val str = """ quote(my_array) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_array")))
    }
  }

  it should "discover an array variable lookup within a squote() call" in {
    val str = """ squote(my_array) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_array")))
    }
  }

  it should "discover an array variable lookup within a unzip() call" in {
    val str = """ unzip(my_array_of_pairs) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("my_array_of_pairs")))
    }
  }

  it should "discover an array variable lookup within a struct literal member access" in {
    val str = """ (StructWithAnArray{myArrayMember: arrayToLookup}).myArray """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.expressionConsumedValueHooks should be(Set(UnlinkedIdentifierHook("arrayToLookup")))
    }
  }
}
