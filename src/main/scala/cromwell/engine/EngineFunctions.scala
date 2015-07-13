package cromwell.engine

import cromwell.binding.WdlFunctions
import cromwell.binding.values._

import scala.util.{Failure, Try}

trait EngineFunctions extends WdlFunctions {
  protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray]
  protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger]
  protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString]
  protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile]
  protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile]

  /**
   * Extract a single `WdlValue` from the specified `Seq`, returning `Failure` if the parameters
   * represent something other than a single `WdlValue`.
   */
  protected def extractSingleArgument(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
    if (params.length != 1) Failure(new UnsupportedOperationException("Expected one argument, got " + params.length))
    else params.head
  }

  def getFunction(name: String): WdlFunction = {
    name match {
      case "read_lines" => read_lines
      case "read_int" => read_int
      case "read_string" => read_string
      case "stdout" => stdout
      case "stderr" => stderr
    }
  }
}
