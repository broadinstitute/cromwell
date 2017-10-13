package cwl

import cats.data.Validated.Valid
import cats.implicits._
import cwl.CommandOutputBinding.Glob
import cwl.ExpressionEvaluator._
import eu.timepit.refined._
import eu.timepit.refined.string.MatchesRegex
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wom.expression.PlaceholderIoFunctionSet
import wom.types.WdlIntegerType
import wom.values.{WdlGlobFile, WdlInteger, WdlString}

class CommandOutputExpressionSpec extends FlatSpec with Matchers {

  behavior of "CommandOutputExpression"

  it should "evaluateValue" in {
    val tempFile = better.files.File.newTemporaryFile("glob.", ".txt").write("41.1")
    val globExpression = Coproduct[Expression](refineMV[MatchesRegex[ECMAScriptExpressionWitness.T]]("$(inputs.myTempFile)"))
    val outputEvalExpression = Coproduct[Expression](refineMV[MatchesRegex[ECMAScriptExpressionWitness.T]]("$((parseInt(self[0].contents) + 1).toFixed())"))
    val glob = Coproduct[Glob](globExpression)
    val outputEval = Coproduct[StringOrExpression](outputEvalExpression)
    val outputBinding = CommandOutputBinding(Option(glob), Option(true), Option(outputEval))
    val commandOutputExpression = CommandOutputExpression(outputBinding, WdlIntegerType)
    val inputValues = Map("myTempFile" -> WdlString(tempFile.pathAsString))
    val result = commandOutputExpression.evaluateValue(inputValues, PlaceholderIoFunctionSet)
    result should be(WdlInteger(42).valid)
  }

  it should "figure out stdout" in {
    val glob = Coproduct[Glob]("stdout")
    val outputBinding = CommandOutputBinding(Option(glob))
    val commandOutputExpression = CommandOutputExpression(outputBinding, WdlIntegerType)
    val result = commandOutputExpression.evaluateFiles(Map.empty, PlaceholderIoFunctionSet, WdlIntegerType)
    result shouldBe(Valid(Set(WdlGlobFile("stdout"))))
  }
}
