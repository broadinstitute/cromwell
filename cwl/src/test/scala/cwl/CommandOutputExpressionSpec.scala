package cwl

import cats.implicits._
import cwl.CommandLineTool.CommandOutputParameter
import cwl.ExpressionEvaluator._
import eu.timepit.refined._
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import wom.expression.{IoFunctionSet, PlaceholderIoFunctionSet}
import wom.types.WomIntegerType
import wom.values._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

class CommandOutputExpressionSpec extends FlatSpec with Matchers {

  behavior of "CommandOutputExpression"

  def ioFunctionSet(data: String) =
    new IoFunctionSet {
      override def readFile(path: String, maxBytes: Option[Int] = None, failOnOverflow: Boolean = false) = Future.successful(data)
      override def writeFile(path: String, content: String) = fail("writeFile should not be used in this test")
      override def copyFile(pathFrom: String, pathTo: String): Future[WomSingleFile] =
        throw new Exception("copyFile should not be used in this test")
      override def stdout(params: Seq[Try[WomValue]]) = fail("stdout should not be used in this test")
      override def stderr(params: Seq[Try[WomValue]]) = fail("stderr should not be used in this test")
      // For this test just "match" the submitted file.
      override def glob(pattern: String): Future[Seq[String]] = Future.successful(List(pattern))
      override def listAllFilesUnderDirectory(dirPath: String): Nothing =
        fail("listAllFilesUnderDirectory should not be used in this test")
      override def size(path: String): Future[Long] = fail("size should not be used in this test")
      override implicit def ec: ExecutionContext = fail("ec should not be used in this test")
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
    val commandOutputExpression = OutputParameterExpression(commandOutputParameter, WomIntegerType, Set.empty, Vector.empty)
    val inputValues = Map("myTempFile" -> WomString(tempFile.pathAsString))
    val result = commandOutputExpression.evaluateValue(inputValues, ioFunctionSet(data))
    result should be(WomInteger(42).valid)
  }

  it should "figure out stdout" in {
    val glob = Coproduct[Glob](Coproduct[StringOrExpression]("stdout"))
    val outputBinding = CommandOutputBinding(Option(glob))
    val commandOutputParameter = CommandOutputParameter("id", outputBinding = Option(outputBinding))
    val commandOutputExpression = OutputParameterExpression(commandOutputParameter, WomIntegerType, Set.empty, Vector.empty)
    val result = commandOutputExpression.evaluateFiles(Map.empty, PlaceholderIoFunctionSet, WomIntegerType)
    // TODO: This should be a glob. See [[cwl.CommandOutputBinding.dependingIfGlobLike]]
    result shouldBe Set(WomSingleFile("stdout")).valid
  }
}
