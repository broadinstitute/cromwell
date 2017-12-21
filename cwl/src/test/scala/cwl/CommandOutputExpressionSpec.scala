package cwl

import cats.data.Validated.Valid
import cats.implicits._
import cwl.CommandOutputBinding.Glob
import cwl.ExpressionEvaluator._
import eu.timepit.refined._
import eu.timepit.refined.string.MatchesRegex
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wom.expression.{IoFunctionSet, PlaceholderIoFunctionSet}
import wom.types.WomIntegerType
import wom.values._

import scala.concurrent.Future
import scala.util.Try

class CommandOutputExpressionSpec extends FlatSpec with Matchers {

  behavior of "CommandOutputExpression"

  def ioFunctionSet(data: String) =
    new IoFunctionSet {
      override def readFile(path: String) = Future.successful(data)
      override def writeFile(path: String, content: String) = throw new Exception("writeFile should not be used in this test")
      override def stdout(params: Seq[Try[WomValue]]) = throw new Exception("stdout should not be used in this test")
      override def stderr(params: Seq[Try[WomValue]]) = throw new Exception("stderr should not be used in this test")
      override def glob(pattern: String): Future[Seq[String]] = throw new Exception("glob should not be used in this test")
      override def size(params: Seq[Try[WomValue]]) = throw new Exception("size should not be used in this test")
    }

  it should "evaluateValue" in {
    val data = "41.1"
    val tempFile = better.files.File.newTemporaryFile("glob.", ".txt").write(data)
    val globExpression = Coproduct[Expression](refineMV[MatchesRegex[ECMAScriptExpressionWitness.T]]("$(inputs.myTempFile)"))
    val outputEvalExpression = Coproduct[Expression](refineMV[MatchesRegex[ECMAScriptExpressionWitness.T]]("$((parseInt(self[0].contents) + 1).toFixed())"))
    val glob = Coproduct[Glob](globExpression)
    val outputEval = Coproduct[StringOrExpression](outputEvalExpression)
    val outputBinding = CommandOutputBinding(Option(glob), Option(true), Option(outputEval))
    val commandOutputExpression = CommandOutputExpression(outputBinding, WomIntegerType, Set.empty)
    val inputValues = Map("myTempFile" -> WomString(tempFile.pathAsString))
    val result = commandOutputExpression.evaluateValue(inputValues, ioFunctionSet(data))
    result should be(WomInteger(42).valid)
  }

  it should "figure out stdout" in {
    val glob = Coproduct[Glob]("stdout")
    val outputBinding = CommandOutputBinding(Option(glob))
    val commandOutputExpression = CommandOutputExpression(outputBinding, WomIntegerType, Set.empty)
    val result = commandOutputExpression.evaluateFiles(Map.empty, PlaceholderIoFunctionSet, WomIntegerType)
    result shouldBe Valid(Set(WomGlobFile("stdout")))
  }
}
