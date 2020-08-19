package wdl.draft3.transforms.expression.consumed

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft3.transforms.linking.expression.consumed.expressionElementUnlinkedValueConsumer
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.UnlinkedIdentifierHook
import wom.values.WomInteger


class ValueConsumerSpec extends AnyFlatSpec with Matchers {
  "the glob value consumer" should "find consumed lookup 'x'" in {
    val expr: ExpressionElement = Glob(IdentifierLookup("x"))
    expr.expressionConsumedValueHooks shouldBe Set(UnlinkedIdentifierHook("x"))
  }

  "the read_string value consumer" should "find consumed lookup 'x' in read_string(x[0])" in {
    val expr: ExpressionElement = ReadString(IndexAccess(IdentifierLookup("x"), PrimitiveLiteralExpressionElement(WomInteger(0))))
    expr.expressionConsumedValueHooks shouldBe Set(UnlinkedIdentifierHook("x"))
  }
}
