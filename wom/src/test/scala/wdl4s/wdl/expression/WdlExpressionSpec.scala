package wdl4s.wdl.expression

import wdl4s.wdl.WdlExpression
import org.scalatest.{FlatSpec, Matchers}

class WdlExpressionSpec extends FlatSpec with Matchers {
  val expr: String => WdlExpression = WdlExpression.fromString

  /* String-ification */
  "Expression Evaluator string-ifier" should "Make strings out of + expressions" in {
    expr("1 + 2").toWdlString shouldEqual "1 + 2"
  }
  it should "Make strings out of - expressions" in {
    expr("1 - 2").toWdlString shouldEqual "1 - 2"
  }
  it should "Make strings out of * expressions" in {
    expr("1 * 2").toWdlString shouldEqual "1 * 2"
  }
  it should "Make strings out of / expressions" in {
    expr("1 / 2").toWdlString shouldEqual "1 / 2"
  }
  it should "Make strings out of % expressions" in {
    expr("1 % 2").toWdlString shouldEqual "1 % 2"
  }
  it should "Make strings out of < expressions" in {
    expr("1 < 2").toWdlString shouldEqual "1 < 2"
  }
  it should "Make strings out of <= expressions" in {
    expr("1 <= 2").toWdlString shouldEqual "1 <= 2"
  }
  it should "Make strings out of > expressions" in {
    expr("1 > 2").toWdlString shouldEqual "1 > 2"
  }
  it should "Make strings out of >= expressions" in {
    expr("1 >= 2").toWdlString shouldEqual "1 >= 2"
  }
  it should "Make strings out of == expressions" in {
    expr("1 == 2").toWdlString shouldEqual "1 == 2"
  }
  it should "Make strings out of != expressions" in {
    expr("1 != 2").toWdlString shouldEqual "1 != 2"
  }
  it should "Make strings out of && expressions" in {
    expr("1 && 2").toWdlString shouldEqual "1 && 2"
  }
  it should "Make strings out of || expressions" in {
    expr("1 || 2").toWdlString shouldEqual "1 || 2"
  }
  it should "Make strings out of expression with strings in it" in {
    expr("\"a\" + \"b\"").toWdlString shouldEqual "\"a\" + \"b\""
  }
  it should "Make strings out of expression with floats in it" in {
    expr("1.1 + 2.2").toWdlString shouldEqual "1.1 + 2.2"
  }
  it should "Make strings out of expression with identifiers in it" in {
    expr("foo + bar").toWdlString shouldEqual "foo + bar"
  }
  it should "Make strings out of member access expressions" in {
    expr("a.b.c").toWdlString shouldEqual "a.b.c"
  }
  it should "Make strings out of function calls" in {
    expr("a(b, c)").toWdlString shouldEqual "a(b, c)"
  }
  it should "Make strings out of array/map lookups" in {
    expr("a[0]").toWdlString shouldEqual "a[0]"
  }
  it should "Make strings out of unary minus" in {
    expr("-2").toWdlString shouldEqual "-2"
  }
  it should "Make strings out of unary plus" in {
    expr("+2").toWdlString shouldEqual "+2"
  }
  it should "Make strings out of logical not" in {
    expr("!2").toWdlString shouldEqual "!2"
  }
  it should "Make strings out of booleans" in {
    expr("true   !=   false").toWdlString shouldEqual "true != false"
  }
}
