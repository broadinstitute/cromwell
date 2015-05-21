package cromwell.engine.backend.local

import cromwell.binding.values.{WdlFile, WdlInteger, WdlString, WdlValue}
import cromwell.engine.EngineFunctions
import cromwell.util.FileUtil

import scala.util.{Failure, Try}


class LocalEngineFunctions(executionContext: TaskExecutionContext) extends EngineFunctions {

  /**
   * Extract a single `WdlValue` from the specified `Seq`, returning `Failure` if the parameters
   * represent something other than a single `WdlValue`.
   */
  private def extractSingleArgument(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
    if (params.length != 1) Failure(new UnsupportedOperationException("Expected one argument, got " + params.length))
    else params.head
  }

  /**
   * Read the entire contents of a file from the specified `WdlValue`, where the file can be
   * specified either as a path via a `WdlString` (with magical handling of "stdout"), or
   * directly as a `WdlFile`.
   *
   * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
   *                                       to be wrapped in a `Try`.
   */
  private def fileContentsToString(value: WdlValue): String = {
    value match {
      case f: WdlFile => FileUtil.slurp(f.value)
      case s: WdlString if s.value == "stdout" => FileUtil.slurp(executionContext.stdout.toFile)
      case e => throw new UnsupportedOperationException("Unsupported argument " + e)
    }
  }

  /**
   * Try to read a string from the file referenced by the specified `WdlValue`.
   */
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument)
    } yield WdlString(string)
  }

  /**
   * Try to read an integer from the file referenced by the specified `WdlValue`.
   */
  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] =
    read_string(params).map { s => WdlInteger(s.value.trim.toInt) }
}
