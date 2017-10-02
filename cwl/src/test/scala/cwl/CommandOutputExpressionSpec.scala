package cwl

import cats.implicits._
import cwl.CommandOutputBinding.Glob
import eu.timepit.refined._
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wdl.types.WdlIntegerType
import wdl.values.{WdlInteger, WdlString}
import wom.expression.PlaceholderIoFunctionSet

class CommandOutputExpressionSpec extends FlatSpec with Matchers {

  behavior of "CommandOutputExpression"

  it should "evaluateValue" in {
    val tempFile = better.files.File.newTemporaryFile("glob.", ".txt").write("41.1")
    val globExpression = refineMV[MatchesECMAScript]("$(inputs.myTempFile)")
    val outputEvalExpression = refineMV[MatchesECMAScript]("$((parseInt(self[0].contents) + 1).toFixed())")
    val glob = Coproduct[Glob](globExpression)
    val outputEval = Coproduct[StringOrExpression](outputEvalExpression)
    val outputBinding = CommandOutputBinding(Option(glob), Option(true), Option(outputEval))
    val commandOutputParameter = CommandOutputParameter("fake_id", outputBinding = Option(outputBinding))
    val commandOutputExpression = CommandOutputExpression(commandOutputParameter, WdlIntegerType)
    val inputValues = Map("myTempFile" -> WdlString(tempFile.pathAsString))
    val result = commandOutputExpression.evaluateValue(inputValues, PlaceholderIoFunctionSet)
    result should be(WdlInteger(42).valid)
  }
}
