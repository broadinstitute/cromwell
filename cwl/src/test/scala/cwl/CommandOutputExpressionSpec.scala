package cwl

import java.util.concurrent.Executors

import cats.syntax.all._
import common.assertion.CromwellTimeoutSpec
import cwl.CommandLineTool.CommandOutputParameter
import cwl.ExpressionEvaluator._
import eu.timepit.refined._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import shapeless.Coproduct
import wom.expression.{EmptyIoFunctionSet, FileEvaluation, NoIoFunctionSet}
import wom.types.WomIntegerType
import wom.values._

import scala.concurrent.{ExecutionContext, Future}


class CommandOutputExpressionSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "CommandOutputExpression"

  def ioFunctionSet(data: String) =
    new EmptyIoFunctionSet {
      override implicit def ec: ExecutionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(1))
      override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false) = Future.successful(data)
      // For this test just "match" the submitted file.
      override def glob(pattern: String): Future[Seq[String]] = Future.successful(List(pattern))
    }

  it should "evaluateValue" in {
    val data = "41.1"
    val tempFile = better.files.File.newTemporaryFile("glob.", ".txt").write(data)
    val globExpression = Coproduct[Expression](refineMV[MatchesECMAScriptExpression]("$(inputs.myTempFile)"))
    val outputEvalExpression = Coproduct[Expression](refineMV[MatchesECMAScriptExpression](
      "$((parseInt(self[0].contents) + 1).toFixed())"))
    val glob = Coproduct[Glob](Coproduct[StringOrExpression](globExpression))
    val outputEval = Coproduct[StringOrExpression](outputEvalExpression)
    val outputBinding = CommandOutputBinding(Option(glob), Option(true), Option(outputEval))
    val commandOutputParameter = CommandOutputParameter("id", outputBinding = Option(outputBinding))
    val commandOutputExpression = OutputParameterExpression(commandOutputParameter, WomIntegerType, Set.empty, Vector.empty, SchemaDefRequirement())
    val inputValues = Map("myTempFile" -> WomString(tempFile.pathAsString))
    val result = commandOutputExpression.evaluateValue(inputValues, ioFunctionSet(data))
    result should be(WomInteger(42).valid)
  }

  it should "figure out stdout" in {
    val glob = Coproduct[Glob](Coproduct[StringOrExpression]("stdout"))
    val outputBinding = CommandOutputBinding(Option(glob))
    val commandOutputParameter = CommandOutputParameter("id", outputBinding = Option(outputBinding))
    val commandOutputExpression = OutputParameterExpression(commandOutputParameter, WomIntegerType, Set.empty, Vector.empty, SchemaDefRequirement())
    val result = commandOutputExpression.evaluateFiles(Map.empty, NoIoFunctionSet, WomIntegerType)
    // TODO: This should be a glob. See [[cwl.CommandOutputBinding.dependingIfGlobLike]]
    result shouldBe Set(FileEvaluation.requiredFile(WomSingleFile("stdout"))).valid
  }
}
