package cwl

import common.validation.Validation._
import cwl.ExpressionEvaluator._
import eu.timepit.refined._
import eu.timepit.refined.numeric.Positive
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wom.callable.RuntimeEnvironment
import wom.expression.NoIoFunctionSet
import wom.graph.LocalName
import wom.values.WomString

class CwlExpressionCommandPartSpec extends FlatSpec with Matchers {

  behavior of "CwlExpressionCommandPart"

  val emptyEnvironment = RuntimeEnvironment("","",refineMV[Positive](1),1,1,1)

  it should "instantiate" in {
    // NOTE: toFixed used to remove the fraction part of ECMAScript numbers
    // https://stackoverflow.com/questions/25989642/why-does-java-8-nashorn-javascript-modulo-returns-0-0-double-instead-of-0-i#answer-25991982
    // https://community.apigee.com/questions/33936/javascript-parseint-not-converting-to-int-value-ne.html
    val commandPart = CwlExpressionCommandPart(Coproduct[Expression](refineMV[MatchesECMAScriptExpression](
      "$(parseInt(inputs.myStringInt).toFixed())"
    )))(false, Vector.empty)
    val result = commandPart.instantiate(Map(LocalName("myStringInt") -> WomString("3")), NoIoFunctionSet, identity, emptyEnvironment).toTry.get.head.commandString
    result should be("'3'")
  }

}
