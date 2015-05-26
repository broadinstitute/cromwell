package cromwell.engine

import cromwell.binding.WdlFunctions
import cromwell.binding.values.{WdlInteger, WdlString, WdlValue}

import scala.util.Try

trait EngineFunctions extends WdlFunctions {
  protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger]
  protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString]

  def getFunction(name: String): WdlFunction = {
    name match {
      case "read_int" => read_int
      case "read_string" => read_string
    }
  }
}
