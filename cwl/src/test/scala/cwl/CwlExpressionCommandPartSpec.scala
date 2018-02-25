package cwl

import common.validation.Validation._
import org.scalatest.{FlatSpec, Matchers}
import eu.timepit.refined._
import eu.timepit.refined.string.MatchesRegex
import ExpressionEvaluator._
import shapeless.Coproduct
import wom.callable.RuntimeEnvironment
import wom.expression.PlaceholderIoFunctionSet
import wom.graph.LocalName
import wom.values.WomString

class CwlExpressionCommandPartSpec extends FlatSpec with Matchers {

  behavior of "CwlExpressionCommandPart"

  val emptyEnvironment = RuntimeEnvironment("","",1,1,1,1)

  it should "instantiate" in {
    // NOTE: toFixed used to remove the fraction part of ECMAScript numbers
    // https://stackoverflow.com/questions/25989642/why-does-java-8-nashorn-javascript-modulo-returns-0-0-double-instead-of-0-i#answer-25991982
    // https://community.apigee.com/questions/33936/javascript-parseint-not-converting-to-int-value-ne.html
    val commandPart = CwlExpressionCommandPart(Coproduct[Expression](refineMV[MatchesRegex[ECMAScriptExpressionWitness.T]]("$(parseInt(inputs.myStringInt).toFixed())")))(Vector.empty)
    val result = commandPart.instantiate(Map(LocalName("myStringInt") -> WomString("3")), PlaceholderIoFunctionSet, identity, emptyEnvironment).toTry.get.head.commandString
    result should be("'3'")
  }

}
