package cwl

import common.validation.Validation._
import org.scalatest.{FlatSpec, Matchers}
import wom.expression.DefaultSizeIoFunctionSet
import wom.values.{WomBoolean, WomString, WomValue}

class ExpressionEvaluatorSpec extends FlatSpec with Matchers {

  behavior of "ExpressionEvaluator"

  it should "eval inputs of different types" in {
    val values: Map[String, WomValue] = Map(
      "myName" -> WomString("hi"),
      "someother" -> WomBoolean(false)
    )
    val parameterContext = ParameterContext(DefaultSizeIoFunctionSet, Vector.empty, values)
    ExpressionEvaluator.eval("inputs.myName", parameterContext).toTry.get should be(WomString("hi"))
    ExpressionEvaluator.eval("inputs.someother", parameterContext).toTry.get should be(WomBoolean(false))
  }
}
