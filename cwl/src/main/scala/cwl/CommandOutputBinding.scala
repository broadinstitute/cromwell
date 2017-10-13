package cwl

import cwl.CommandOutputBinding.Glob
import shapeless.{:+:, CNil}
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

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
  def commandOutputBindingToWdlValue(parameterContext: ParameterContext,
                                     ioFunctionSet: IoFunctionSet): WdlValue = {

    val paths: Seq[String] = glob map { globValue =>
      GlobEvaluator.globPaths(globValue, parameterContext, ioFunctionSet)
    } getOrElse {
      Vector.empty
    }

    val loadContents: Boolean = this.loadContents getOrElse false

    val wdlMapType = WdlMapType(WdlStringType, WdlStringType)
    val wdlMaps = paths map { path =>
      // TODO: WOM: basename/dirname/size/checksum/etc.

      val contents: Map[WdlValue, WdlValue] =
        if (loadContents) Map(WdlString("contents") -> WdlString(load64KiB(path, ioFunctionSet))) else Map.empty

      val wdlKeyValues: Map[WdlValue, WdlValue] = Map(
        WdlString("location") -> WdlString(path)
      ) ++ contents

      WdlMap(wdlMapType, wdlKeyValues)
    }

    val arrayOfCwlFileMaps = WdlArray(WdlArrayType(wdlMapType), wdlMaps)

    val outputEvalParameterContext = parameterContext.copy(self = arrayOfCwlFileMaps)

    outputEval match {
      case Some(outputEvalCoproduct) =>
        outputEvalCoproduct.fold(OutputEvalToWdlValue).apply(outputEvalParameterContext)
      case None =>
        // Return the WdlArray of file paths, three_step.ps needs this for stdout output.
        // There will be conversion required between this Array[File] output type and the requested File.
        arrayOfCwlFileMaps
    }
  }


  private def load64KiB(path: String, ioFunctionSet: IoFunctionSet): String = {
    // This suggests the IoFunctionSet should have a length-limited read API as both CWL and WDL support this concept.
    // val content = ioFunctionSet.readFile(path)
    val content = better.files.File(path).bytes.take(64 * 1024).toArray
    new String(content)
  }
}

object CommandOutputBinding {
  type Glob = Expression :+: String :+: Array[String] :+: CNil

}
