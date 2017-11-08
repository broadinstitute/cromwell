package cwl

import cwl.CommandOutputBinding.Glob
import shapeless.{:+:, CNil}
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

import scala.language.postfixOps
import scala.concurrent.Await
import scala.concurrent.duration._

/** @see <a href="http://www.commonwl.org/v1.0/Workflow.html#CommandOutputBinding">CommandOutputBinding</a> */
case class CommandOutputBinding(
                                 glob: Option[Glob] = None,
                                 loadContents: Option[Boolean] = None,
                                 outputEval: Option[StringOrExpression] = None) {

  /*
  CommandOutputBinding.glob:
  Find files relative to the output directory, using POSIX glob(3) pathname matching. If an array is provided, find
  files that match any pattern in the array. If an expression is provided, the expression must return a string or an
  array of strings, which will then be evaluated as one or more glob patterns. Must only match and return files which
  actually exist.

  http://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputBinding
   */
  def commandOutputBindingToWomValue(parameterContext: ParameterContext,
                                     ioFunctionSet: IoFunctionSet): WomValue = {

    val paths: Seq[String] = glob map { globValue =>
      GlobEvaluator.globPaths(globValue, parameterContext, ioFunctionSet)
    } getOrElse {
      Vector.empty
    }

    val loadContents: Boolean = this.loadContents getOrElse false

    val wdlMapType = WomMapType(WomStringType, WomStringType)
    val wdlMaps = paths map { path =>
      // TODO: WOM: basename/dirname/size/checksum/etc.

      val contents: Map[WomValue, WomValue] =
        if (loadContents) Map(WomString("contents") -> WomString(load64KiB(path, ioFunctionSet))) else Map.empty

      val wdlKeyValues: Map[WomValue, WomValue] = Map(
        WomString("location") -> WomString(path)
      ) ++ contents

      WomMap(wdlMapType, wdlKeyValues)
    }

    val arrayOfCwlFileMaps = WomArray(WomArrayType(wdlMapType), wdlMaps)

    val outputEvalParameterContext = parameterContext.copy(self = arrayOfCwlFileMaps)

    outputEval match {
      case Some(outputEvalCoproduct) =>
        outputEvalCoproduct match {
          case StringOrExpression.String(s) => WomString(s)
          case StringOrExpression.Expression(e) => e.fold(EvaluateExpression).apply(outputEvalParameterContext)
        }
      case None =>
        // Return the WdlArray of file paths, three_step.ps needs this for stdout output.
        // There will be conversion required between this Array[File] output type and the requested File.
        arrayOfCwlFileMaps
    }
  }


  private def load64KiB(path: String, ioFunctionSet: IoFunctionSet): String = {
    // This suggests the IoFunctionSet should have a length-limited read API as both CWL and WDL support this concept.
    // ChrisL: But remember that they are different (WDL => failure, CWL => truncate)
    val content = ioFunctionSet.readFile(path)

    // TODO: propagate IO, Try, or Future or something all the way out via "commandOutputBindingtoWdlValue" signature
    // TODO: Stream only the first 64 KiB, this "read everything then ignore most of it" method is terrible
    val initialResult = Await.result(content, 5 seconds)
    initialResult.substring(0, Math.min(initialResult.length, 64 * 1024))
  }
}

object CommandOutputBinding {
  type Glob = Expression :+: String :+: Array[String] :+: CNil

}
